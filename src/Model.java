import org.apache.commons.math3.analysis.interpolation.*;

/**
 * Created by lahmann on 2016-06-19.
 *
 * Class the holds all of the information required to sample particles from some ParticleDistribution
 * and track them through some Plasma. The actual step by step logic of the simulation is handled within
 * individual ParticleHistoryTasks, this class just initiates them
 *
 */
public class Model {

    private ParticleDistribution testParticleDistribution;
    private StoppingPowerModel stoppingPower;
    private CrossSection crossSection;
    private Plasma plasma;



    public Model(ParticleDistribution testParticleDistribution, CrossSection crossSection, Plasma plasma) {
        this.testParticleDistribution = testParticleDistribution;
        this.crossSection = crossSection;
        this.plasma = plasma;
        this.stoppingPower = new StoppingPowerModel(testParticleDistribution, plasma);
    }

    public double getYieldRatio(int totalParticles){

        int numCPUs = Runtime.getRuntime().availableProcessors();
        //int numCPUs = 1;
        int totalPerThread = (int) Math.ceil((double) totalParticles / numCPUs);

        Thread[] threads = new Thread[numCPUs];
        ParticleHistoryTask[] tasks = new ParticleHistoryTask[numCPUs];


        for (int i = 0; i < numCPUs; i++){
            tasks[i] = new ParticleHistoryTask(
                    crossSection, testParticleDistribution, plasma, stoppingPower, totalPerThread);

            //tasks[i].setDebugMode(true);
            //tasks[i].run();

            threads[i] = new Thread(tasks[i]);
            threads[i].start();
        }

        double yieldRatio = 0.0;
        for (int i = 0; i < numCPUs; i ++){
            try {
                threads[i].join();
                yieldRatio += tasks[i].getReactionProbability();

            }catch (Exception e){
                e.printStackTrace();
            }
        }

        return yieldRatio/numCPUs;
    }

}
