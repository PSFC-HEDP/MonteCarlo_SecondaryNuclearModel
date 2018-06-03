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

    // Reactivity of the source particle (used to establish the spatial profile)
    private Reactivity sourceReactivity;

    // An ORDERED array list of all the plasmas in this model
    private ArrayList<Plasma> plasmas = new ArrayList<>();

    // An array list of all nuclear reactions we'll simulate in this model
    private ArrayList<NuclearReaction> nuclearReactions = new ArrayList<>();

    // The line of sight where tallies will be done (can be null) TODO: Consider an array list for multiple tallys
    private Vector3D detectorLineOfSight;

    // Threads that this model interacts with
    private SimulationThread[] threads;

    // Logger object that will keep track of progress / errors
    private Logger logger = new Logger(true);

    // Array of tallies for the source particles (1 tally for birth, and 1 for every surface)
    private Tally[] sourceParticlePositionTallies = null;
    private Tally[] sourceParticleEnergyTallies = null;
    private Tally[] sourceParticleTimeTallies = null;

    // HashMap of tally arrays for every nuclear reaction we model
    private HashMap<NuclearReaction, Tally[]> productParticlePositionTallyMap = null;
    private HashMap<NuclearReaction, Tally[]> productParticleEnergyTallyMap = null;
    private HashMap<NuclearReaction, Tally[]> productParticleTimeTallyMap = null;



    /**
     * Default constructor
     * @param name Filename for the file where the output data will be saved
     */
    public Model(String name) {
        this.outputFile = new OutputFile(name);
    }

    /**
     * Preconfigured model to simulate DD-tritions and generate secondary DT neutrons
     * User still must specify a source plasma with deuterium
     * @param name Filename for the file where the output data will be saved
     * @return A Model for secondary DTn simulations
     */
    public static Model SecondaryDTnModel(String name){
        Model model = new Model(name);
        model.setSourceInformation(NuclearReaction.DD_t, Reactivity.DDp_Reactivity);
        model.addNuclearReaction(NuclearReaction.DT_n);
        return model;
    }


    /**
     * Method to set all of the required source information. This method much be called to have a valid model
     * @param sourceNuclearReaction Nuclear reaction of the source (primary) particle
     * @param sourceReactivity Reactivity of the source particle reaction
     */
    public void setSourceInformation(NuclearReaction sourceNuclearReaction, Reactivity sourceReactivity) {
        this.sourceNuclearReaction = sourceNuclearReaction;
        this.sourceReactivity = sourceReactivity;
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
     * Method to set the physical location of the tallys. If not set, tallys will be averaged over 4pi.
     * @param detectorLineOfSight Physical tally location
     */
    public void setDetectorLineOfSight(Vector3D detectorLineOfSight) {
        this.detectorLineOfSight = detectorLineOfSight;
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

            // Verify that the source information has been set
            if (sourceNuclearReaction == null || sourceReactivity == null)
                throw new Exceptions.NoSourceInformationException();


            // Verify that we have a plasma
            if (plasmas.size() == 0)
                throw new Exceptions.NoPlasmaSpecifiedException();


            // Verify this machine has enough processors
            if (numCPUs > Runtime.getRuntime().availableProcessors() - 1 || numCPUs < 1)
                throw new Exceptions.InvalidNumberProcessorsRequestedException(numCPUs, Runtime.getRuntime().availableProcessors() - 1);


            // Calculate the number of particles we'll simulate in each thread
            int particlesPerThread = Math.floorDiv(totalParticles, numCPUs) + 1;


            // Organize the all of the plasma data to save runtime
            logger.addLog("Starting model setup ... ");
            PlasmaData[] plasmaData = new PlasmaData[plasmas.size()];
            for (int i = 0; i < plasmaData.length; i++) {
                plasmaData[i] = new PlasmaData(plasmas.get(i));
            }
            logger.logTaskCompletion();



            /**
             * Running the actual simulation
             */

            // Set up the simulation threads
            logger.addLog(String.format("Building %d threads with %d particles each ... ", numCPUs, particlesPerThread));
            threads = new SimulationThread[numCPUs];
            for (int i = 0; i < threads.length; i++) {

                // Make the thread
                threads[i] = new SimulationThread(particlesPerThread);
                logger.addLog("Started thread " + threads[i].getId());

                // Add the source and tally info
                threads[i].setSourceInformation(plasmas.get(0), sourceReactivity, sourceNuclearReaction);
                threads[i].setDetectorLineOfSight(detectorLineOfSight);

                // Add the plasma data
                for (PlasmaData data : plasmaData) {
                    threads[i].addPlasma(data);
                }

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


        /*
        // Start the output file
        FileWriter w = new FileWriter(outputFile);

        // Write the plasma information
        for (int i = 0; i < plasmas.size(); i++){
            w.write("Plasma Layer " + i + "\n");
            w.write(plasmas.get(i).toString() + "\n");
        }

        // Write the source distribution to the output file
        w.write("Source Particle Nuclear Reaction" + "\n");
        w.write(sourceNuclearReaction.toString() + "\n");

        // Write all of the modeled nuclear reactions to the output file
        for (int i = 0; i < nuclearReactions.size(); i++){
            w.write("Modelled Nuclear Reaction " + i + "\n");
            w.write(nuclearReactions.get(i).toString() + "\n");
        }

        w.write("Job started at " + System.currentTimeMillis() + "\n");


        // Set up the threads and tasks
        Thread[] threads = new Thread[numCPUs];
        SimulationThread[] tasks = new SimulationThread[numCPUs];
        for (int i = 0; i < numCPUs; i++){

            // Sanity check to make sure we don't simulate too many and mess up our normalization
            int particlesForTask = Math.min(totalPerThread, totalParticles);
            totalParticles -= particlesForTask;


            // Create the task
            tasks[i] = new SimulationThread(particlesForTask);
            tasks[i].setSourceInformation(plasmas.get(0), sourceReactivity, sourceNuclearReaction);     // TODO: Assumes the source plasma
            tasks[i].setDetectorLineOfSight(detectorLineOfSight);


            // Sort out all of the data Objects on initiation so we don't have to do it during runtime
            for (Plasma layer : plasmas){
                PlasmaData data = buildLayerData(layer);
                tasks[i].addPlasma(data);
            }


            // Easier for debugging if we run the task in this thread
            //tasks[i].setDebugMode(true);
            //tasks[i].run();


            // Make a note in the output file
            w.write(String.format("%d: Starting task %d with %d particles...\n",
                    System.currentTimeMillis(), i, particlesForTask));


            // Start the task
            threads[i] = new Thread(tasks[i]);
            threads[i].start();
        }


        Tally[] sourceParticlePositionTallies = null;
        Tally[] sourceParticleEnergyTallies   = null;
        Tally[] sourceParticleTimeTallies     = null;

        HashMap<NuclearReaction, Tally[]> productParticlePositionTallyMap = null;
        HashMap<NuclearReaction, Tally[]> productParticleEnergyTallyMap   = null;
        HashMap<NuclearReaction, Tally[]> productParticleTimeTallyMap     = null;

        for (int i = 0; i < numCPUs; i ++){
            try {
                // Join with the thread
                threads[i].join();

                // Grab the source particle tallies
                Tally[] taskSourceParticlePositionTallies = tasks[i].getSourceParticlePositionTallies();
                Tally[] taskSourceParticleEnergyTallies   = tasks[i].getSourceParticleEnergyTallies();
                Tally[] taskSourceParticleTimeTallies     = tasks[i].getSourceParticleTimeTallies();

                // Grab the product particle tally maps
                HashMap<NuclearReaction, Tally[]> taskProductParticlePositionTallyMap = tasks[i].getProductParticlePositionTallyMap();
                HashMap<NuclearReaction, Tally[]> taskProductParticleEnergyTallyMap   = tasks[i].getProductParticleEnergyTallyMap();
                HashMap<NuclearReaction, Tally[]> taskProductParticleTimeTallyMap     = tasks[i].getProductParticleTimeTallyMap();

                // Add the source particle tallies together
                if (i == 0){
                    sourceParticlePositionTallies = taskSourceParticlePositionTallies;
                    sourceParticleEnergyTallies   = taskSourceParticleEnergyTallies;
                    sourceParticleTimeTallies     = taskSourceParticleTimeTallies;
                }else{
                    for (int j = 0; j < sourceParticlePositionTallies.length; j++){
                        sourceParticlePositionTallies[j].addTally(taskSourceParticlePositionTallies[j]);
                        sourceParticleEnergyTallies  [j].addTally(taskSourceParticleEnergyTallies  [j]);
                        sourceParticleTimeTallies    [j].addTally(taskSourceParticleTimeTallies    [j]);
                    }
                }

                // Add the product particle tallies together
                if (i == 0){
                    productParticlePositionTallyMap = taskProductParticlePositionTallyMap;
                    productParticleEnergyTallyMap   = taskProductParticleEnergyTallyMap;
                    productParticleTimeTallyMap     = taskProductParticleTimeTallyMap;
                }else{
                    for (NuclearReaction key : taskProductParticlePositionTallyMap.keySet()){
                        Tally[] taskProductParticlePositionTallies = taskProductParticlePositionTallyMap.get(key);
                        Tally[] taskProductParticleEnergyTallies   = taskProductParticleEnergyTallyMap  .get(key);
                        Tally[] taskProductParticleTimeTallies     = taskProductParticleTimeTallyMap    .get(key);

                        Tally[] productParticlePositionTallies = productParticlePositionTallyMap.get(key);
                        Tally[] productParticleEnergyTallies   = productParticleEnergyTallyMap  .get(key);
                        Tally[] productParticleTimeTallies     = productParticleTimeTallyMap    .get(key);

                        for (int j = 0; j < productParticlePositionTallies.length; j++){
                            productParticlePositionTallies[j].addTally(taskProductParticlePositionTallies[j]);
                            productParticleEnergyTallies  [j].addTally(taskProductParticleEnergyTallies  [j]);
                            productParticleTimeTallies    [j].addTally(taskProductParticleTimeTallies    [j]);
                        }

                        productParticlePositionTallyMap.put(key, productParticlePositionTallies);
                        productParticleEnergyTallyMap  .put(key, productParticleEnergyTallies);
                        productParticleTimeTallyMap    .put(key, productParticleTimeTallies);
                    }
                }

            }catch (Exception e){
                e.printStackTrace();
            }
        }

        w.write(System.currentTimeMillis() + ": All tasks completed!\n\n");

        w.write("Source Particle Tallies\n");
        for (int i = 0; i < sourceParticlePositionTallies.length; i++){

            if (i == 0){
                w.write("Birth Position Tally\n");
                w.write(sourceParticlePositionTallies[i] + "\n");
                w.write("Birth Energy Tally\n");
                w.write(sourceParticleEnergyTallies[i] + "\n");
                w.write("Birth Time Tally\n");
                w.write(sourceParticleTimeTallies[i] + "\n");
            }else{
                w.write("Position Tally Crossing Plasma Layer " + (i-1) + "\n");
                w.write(sourceParticlePositionTallies[i] + "\n");
                w.write("Energy Tally Crossing Plasma Layer " + (i-1) + "\n");
                w.write(sourceParticleEnergyTallies[i] + "\n");
                w.write("Time Tally Crossing Plasma Layer " + (i-1) + "\n");
                w.write(sourceParticleTimeTallies[i] + "\n");

                //double[] fit = sourceParticleTallies[i].getGaussFit();
                //System.out.print(fit[1] + " " + fit[2] + " ");
            }
        }

        w.write("Product Particle Tallies\n");
        for (NuclearReaction key : productParticlePositionTallyMap.keySet()) {
            w.write(key + " Product Particle Tallies\n");

            Tally[] productParticlePositionTallies = productParticlePositionTallyMap.get(key);
            Tally[] productParticleEnergyTallies   = productParticleEnergyTallyMap  .get(key);
            Tally[] productParticleTimeTallies     = productParticleTimeTallyMap    .get(key);

            for (int i = 0; i < productParticlePositionTallies.length; i++) {

                // TODO: TEMP LAZY
                System.out.printf(" %.4e", productParticleEnergyTallies[i].getTotalWeight());

                if (i == 0) {
                    w.write("Birth Position Tally\n");
                    w.write(productParticlePositionTallies[i] + "\n");
                    w.write("Birth Energy Tally\n");
                    w.write(productParticleEnergyTallies[i] + "\n");
                    w.write("Birth Time Tally\n");
                    w.write(productParticleTimeTallies[i] + "\n");
                } else {
                    w.write("Position Tally Crossing Plasma Layer " + (i - 1) + "\n");
                    w.write(productParticlePositionTallies[i] + "\n");
                    w.write("Energy Tally Crossing Plasma Layer " + (i - 1) + "\n");
                    w.write(productParticleEnergyTallies[i] + "\n");
                    w.write("Time Tally Crossing Plasma Layer " + (i - 1) + "\n");
                    w.write(productParticleTimeTallies[i] + "\n");
                }
            }
        }

        System.out.println();
        w.close();
        */
    }

    private void writeOutputFile() throws Exception {

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

    class PlasmaData {

        Plasma plasma;
        HashMap<ParticleType, StoppingPowerModel> stoppingPowerModels = new HashMap<>();
        HashMap<ParticleType, ArrayList<ReactionData>> reactionDataMap = new HashMap<>();


        public PlasmaData(Plasma plasma) {

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
                    reactionDataMap.get(sourceParticle).add(new ReactionData(reaction, B));
                }
            }
        }
    }

    class ReactionData{

        NuclearReaction nuclearReaction;
        ParticleType backgroundParticle;

        public ReactionData(NuclearReaction nuclearReaction, ParticleType backgroundParticle) {
            this.nuclearReaction = nuclearReaction;
            this.backgroundParticle = backgroundParticle;
        }
    }

}
