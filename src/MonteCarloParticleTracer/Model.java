package MonteCarloParticleTracer;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by lahmann on 2016-06-19.
 *
 * Class the holds all of the information required to sample particles from some ParticleDistribution
 * and track them through some PlasmaLayer. The actual step by step logic of the simulation is handled within
 * individual ParticleHistoryTasks, this class just initiates them
 *
 */
public class Model {

    private String name;

    private ParticleDistribution testParticleDistribution;

    private ArrayList<PlasmaLayer> plasmaLayers = new ArrayList<>();
    private ArrayList<NuclearReaction> nuclearReactions = new ArrayList<>();

    private Vector3D detectorLineOfSight;

    public Model(String name) {
        this.name = name;
    }

    public void addPlasmaLayer(PlasmaLayer plasmaLayer){

        // Force the new layer to share a boundary with the out most layer
        if (plasmaLayers.size() > 1){
            PlasmaLayer outerMostLayer = plasmaLayers.get(plasmaLayers.size() - 1);
            plasmaLayer.setInnerBoundaryLegendreModes(outerMostLayer.getOuterBoundaryLegendreModes());
        }

        plasmaLayers.add(plasmaLayer);
    }

    public void addNuclearReaction(NuclearReaction reaction){
        this.addNuclearReaction(reaction);
    }

    public void setTestParticleDistribution(ParticleDistribution testParticleDistribution) {
        this.testParticleDistribution = testParticleDistribution;
    }

    public void setDetectorLineOfSight(Vector3D detectorLineOfSight) {
        this.detectorLineOfSight = detectorLineOfSight;
    }

    public void runSimulation(int totalParticles, int numCPUs, double[] energyNodes){

        // Normalize the distribution
        testParticleDistribution.setParticleWeight(1.0 / totalParticles);

        // Sanity check the number of requested nodes
        if (numCPUs > Runtime.getRuntime().availableProcessors()){
            System.err.printf("Warning! The user requested %d nodes when only %d are available.\n" +
                    "Changing requested number of nodes to %d.\n",
                    numCPUs, Runtime.getRuntime().availableProcessors(), Runtime.getRuntime().availableProcessors());
            numCPUs = Runtime.getRuntime().availableProcessors();
        }

        // Calculate the number of particles per thread
        int totalPerThread = (int) Math.ceil((double) totalParticles / numCPUs);


        // Set up the threads and tasks
        Thread[] threads = new Thread[numCPUs];
        ParticleHistoryTask[] tasks = new ParticleHistoryTask[numCPUs];
        for (int i = 0; i < numCPUs; i++){

            // Sanity check to make sure we don't simulate too many and mess up our normalization
            int particlesForTask = Math.min(totalPerThread, totalParticles);
            totalParticles -= particlesForTask;

            // Create the task
            System.out.printf("Starting task %d with %d particles...\n", i, particlesForTask);
            tasks[i] = new ParticleHistoryTask(testParticleDistribution, particlesForTask);

            // Sort out all of the data Objects on initiation so we don't have to do it during runtime
            for (PlasmaLayer layer : plasmaLayers){
                LayerData data = buildLayerData(testParticleDistribution, layer);
                tasks[i].addPlasmaLayer(data);
            }

            // Easier for debugging if we run the task in this thread
            //tasks[i].setDebugMode(true);
            //tasks[i].run();

            threads[i] = new Thread(tasks[i]);
            threads[i].start();
        }

        // Join with all the threads and add up all the tallies
        Tally escapingParticleEnergyTally = null;
        for (int i = 0; i < numCPUs; i ++){
            try {
                threads[i].join();

                if (i == 0){
                    escapingParticleEnergyTally = tasks[i].getEscapingParticleEnergyTally();
                }else{
                    escapingParticleEnergyTally.addTally(tasks[i].getEscapingParticleEnergyTally());
                }

            }catch (Exception e){
                e.printStackTrace();
            }
        }

        System.out.println(escapingParticleEnergyTally);
    }

    /**
     * Private convenience functions
     */
    private LayerData buildLayerData(ParticleDistribution distribution, PlasmaLayer plasmaLayer){
        LayerData layerData = new LayerData(plasmaLayer, new StoppingPowerModel(distribution, plasmaLayer));

        for (int i = 0; i < nuclearReactions.size(); i++){
            NuclearReaction reaction = nuclearReactions.get(i);

            if (distribution.getType().equals(reaction.getReactantParticleTypeA())){
                if (plasmaLayer.containsSpecies(reaction.getReactantParticleTypeB())){
                    ReactionData data = new ReactionData(reaction, reaction.getReactantParticleTypeB(), i);
                    layerData.addReactionData(data);
                }
            }
            else if (distribution.getType().equals(reaction.getReactantParticleTypeB())){
                if (plasmaLayer.containsSpecies(reaction.getReactantParticleTypeA())){
                    ReactionData data = new ReactionData(reaction, reaction.getReactantParticleTypeA(), i);
                    layerData.addReactionData(data);
                }
            }
        }
        return layerData;
    }

    public class LayerData{

        PlasmaLayer layer;
        StoppingPowerModel stoppingPowerModel;
        ArrayList<ReactionData> reactionData = new ArrayList<>();

        public LayerData(PlasmaLayer layer, StoppingPowerModel stoppingPowerModel) {
            this.layer = layer;
            this.stoppingPowerModel = stoppingPowerModel;
        }

        public void addReactionData(ReactionData data){
            this.reactionData.add(data);
        }
    }

    public class ReactionData{

        NuclearReaction nuclearReaction;
        ParticleType reactingSpecies;
        Integer tallyIndex;

        public ReactionData(NuclearReaction nuclearReaction, ParticleType reactingSpecies, Integer tallyIndex) {
            this.nuclearReaction = nuclearReaction;
            this.reactingSpecies = reactingSpecies;
            this.tallyIndex = tallyIndex;
        }
    }

}
