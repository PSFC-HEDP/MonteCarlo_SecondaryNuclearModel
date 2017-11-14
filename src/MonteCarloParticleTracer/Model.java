package MonteCarloParticleTracer;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.io.File;
import java.io.FileWriter;
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

    private File outputFile;

    private NuclearReaction sourceNuclearReaction;
    private Reactivity sourceReactivity;

    private ArrayList<Plasma> plasmas = new ArrayList<>();
    private ArrayList<NuclearReaction> nuclearReactions = new ArrayList<>();

    private Vector3D detectorLineOfSight;

    public Model(File outputFile) {
        this.outputFile = outputFile;
    }

    public void setSourceInformation(NuclearReaction sourceNuclearReaction, Reactivity sourceReactivity) {
        this.sourceNuclearReaction = sourceNuclearReaction;
        this.sourceReactivity = sourceReactivity;
    }

    public void addPlasmaLayer(Plasma plasma){

        // Force the new layer to share a boundary with the out most layer
        if (plasmas.size() > 0){
            Plasma outerMostLayer = plasmas.get(plasmas.size() - 1);
            plasma.setInnerBoundaryLegendreModes(outerMostLayer.getOuterBoundaryLegendreModes());
        }

        plasmas.add(plasma);
    }

    public void addNuclearReaction(NuclearReaction reaction){
        nuclearReactions.add(reaction);
    }

    public void setDetectorLineOfSight(Vector3D detectorLineOfSight) {
        this.detectorLineOfSight = detectorLineOfSight;
    }

    public void runSimulation(int totalParticles, int numCPUs) throws Exception{

        // Sanity check the number of requested CPUs
        if (numCPUs > Runtime.getRuntime().availableProcessors()){
            System.err.printf("Warning! The user requested %d nodes when only %d are available on this machine.\n" +
                            "Changing requested number of nodes to %d.\n",
                    numCPUs, Runtime.getRuntime().availableProcessors(), Runtime.getRuntime().availableProcessors());
            numCPUs = Runtime.getRuntime().availableProcessors();
        }

        // Calculate the number of particles per thread
        int totalPerThread = (int) Math.ceil((double) totalParticles / numCPUs);


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
        ParticleHistoryTask[] tasks = new ParticleHistoryTask[numCPUs];
        for (int i = 0; i < numCPUs; i++){

            // Sanity check to make sure we don't simulate too many and mess up our normalization
            int particlesForTask = Math.min(totalPerThread, totalParticles);
            totalParticles -= particlesForTask;


            // Create the task
            tasks[i] = new ParticleHistoryTask(particlesForTask);
            tasks[i].setSourceInformation(plasmas.get(0), sourceReactivity, sourceNuclearReaction);     // TODO: Assumes the source plasma
            tasks[i].setDetectorLineOfSight(detectorLineOfSight);


            // Sort out all of the data Objects on initiation so we don't have to do it during runtime
            for (Plasma layer : plasmas){
                LayerData data = buildLayerData(layer);
                tasks[i].addPlasmaLayer(data);
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
    }

    /**
     * Private convenience functions
     */
    private LayerData buildLayerData(Plasma plasma){

        LayerData layerData = new LayerData(plasma);
        ParticleType sourceParticle = sourceNuclearReaction.getProducts()[0];

        // Build the stopping power model for the source particles
        layerData.addStoppingPowerModel(sourceParticle);

        // Loop through all of the nuclear reactions
        for (int i = 0; i < nuclearReactions.size(); i++){

            NuclearReaction reaction = nuclearReactions.get(i);
            ParticleType productParticle = reaction.getProducts()[0];

            // Build the stopping power model for the product of this reaction
            layerData.addStoppingPowerModel(productParticle);

            // Determine which of the reactions can be set off by the source particles
            ParticleType A = reaction.getReactants()[0];
            ParticleType B = reaction.getReactants()[1];

            if (sourceParticle.equals(A) && plasma.containsSpecies(B)){
                layerData.addReactionData(reaction, A, B);
            }

            else if (sourceParticle.equals(B) && plasma.containsSpecies(A)){
                layerData.addReactionData(reaction, B, A);
            }

            /**
             * This is where we would add a loop to check for tertiary reactions if we're ever interested in that
             * We need to consider this carefully though, it may cause severe run time penalties if we don't implement
             * some kind of weight threshold
             */

        }
        return layerData;
    }

    public class LayerData{

        Plasma layer;
        HashMap<ParticleType, StoppingPowerModel> stoppingPowerModels = new HashMap<>();
        HashMap<ParticleType, ArrayList<ReactionData>> reactionDataMap = new HashMap<>();

        public LayerData(Plasma layer) {
            this.layer = layer;
        }

        public void addStoppingPowerModel(ParticleType particleType){
            StoppingPowerModel model = new StoppingPowerModel(particleType, this.layer);
            stoppingPowerModels.put(particleType, model);

            // We'll build an empty entry just in case no reactions get tied to this particle type
            reactionDataMap.put(particleType, new ArrayList<>());
        }

        public void addReactionData(NuclearReaction reaction, ParticleType simulatedParticle, ParticleType backgroundParticle){
            ReactionData data = new ReactionData(reaction, backgroundParticle);

            ArrayList<ReactionData> arrayList = reactionDataMap.get(simulatedParticle);
            // TODO: Potential for double counting here...
            arrayList.add(data);
            reactionDataMap.put(simulatedParticle, arrayList);
        }
    }

    public class ReactionData{

        NuclearReaction nuclearReaction;
        ParticleType backgroundParticle;

        public ReactionData(NuclearReaction nuclearReaction, ParticleType backgroundParticle) {
            this.nuclearReaction = nuclearReaction;
            this.backgroundParticle = backgroundParticle;
        }
    }

}
