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
    private NuclearReaction nuclearReaction;

    private Plasma plasma;
    private ArrayList<StoppingPowerModel> stoppingPowerModels;

    // TODO: TEMP!
    public double[] energyNodes;
    public double[] energyTallies;


    public Model(String name, ParticleDistribution testParticleDistribution, NuclearReaction nuclearReaction, Plasma plasma) {

        this.name = name;

        this.testParticleDistribution = testParticleDistribution;
        this.nuclearReaction = nuclearReaction;

        this.plasma = plasma;

        stoppingPowerModels = new ArrayList<>();
        for (PlasmaLayer layer : plasma.getLayers()){
            stoppingPowerModels.add(new StoppingPowerModel(testParticleDistribution, layer));
        }
    }

    public void runSimulation(int totalParticles, Vector3D detectorLineOfSight, double[] energyNodes){

        // See how many CPUs are available to use (ideally the user should specify this in the future)
        int numCPUs = Runtime.getRuntime().availableProcessors();
        //int numCPUs = 1;
        int totalPerThread = (int) Math.ceil((double) totalParticles / numCPUs);


        // We're rounding up the number of particles to keep the distribution even, so we need to recalculate it
        totalParticles = numCPUs * totalPerThread;
        testParticleDistribution.setParticleWeight(1.0 / totalParticles);


        // Set up the threads and tasks
        Thread[] threads = new Thread[numCPUs];
        ParticleHistoryTask[] tasks = new ParticleHistoryTask[numCPUs];
        for (int i = 0; i < numCPUs; i++){
            tasks[i] = new ParticleHistoryTask(testParticleDistribution, nuclearReaction, plasma.getLayers(),
                    stoppingPowerModels, totalPerThread, detectorLineOfSight, energyNodes);

            //tasks[i].setDebugMode(true);
            //tasks[i].run();


            threads[i] = new Thread(tasks[i]);
            threads[i].start();
        }

        double reactionProbability = 0.0;
        this.energyNodes = energyNodes;
        this.energyTallies = new double[energyNodes.length];

        for (int i = 0; i < numCPUs; i ++){
            try {
                threads[i].join();

                reactionProbability += tasks[i].getReactionProbability();
                ParticleHistoryTask.EnergyTally tally = tasks[i].getEnergyTally();

                for (int j = 0; j < tally.tallies.length; j++){
                    energyTallies[j] += tally.tallies[j];
                }

            }catch (Exception e){
                e.printStackTrace();
            }
        }


        try {
            FileWriter w = new FileWriter(new File(name + "_" + System.currentTimeMillis() + ".csv"));

            w.write(name + "\n");
            w.write("\n");
            w.write(String.format("Total reaction probability:, %.4e\n", reactionProbability));
            w.write("\n");
            w.write("Resultant Spectrum\n");
            w.write("   Nodes   ,   Value\n");
            for (int i = 0; i < energyNodes.length; i++){
                w.write(String.format("%.4e , %.4e\n", energyNodes[i], energyTallies[i]));
            }

            w.close();
        }
        catch (IOException e){
            e.printStackTrace();
        }


    }

}
