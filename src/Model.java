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
    private NuclearReaction nuclearReaction;
    private Plasma plasma;



    public Model(ParticleDistribution testParticleDistribution, NuclearReaction nuclearReaction, Plasma plasma) {
        this.testParticleDistribution = testParticleDistribution;
        this.nuclearReaction = nuclearReaction;
        this.plasma = plasma;
        this.stoppingPower = new StoppingPowerModel(testParticleDistribution, plasma);
    }

    public double getYieldRatio(int totalParticles){

        testParticleDistribution.setParticleWeight(1.0 / totalParticles);

        //int numCPUs = Runtime.getRuntime().availableProcessors();
        int numCPUs = 1;
        int totalPerThread = (int) Math.ceil((double) totalParticles / numCPUs);

        Thread[] threads = new Thread[numCPUs];
        ParticleHistoryTask[] tasks = new ParticleHistoryTask[numCPUs];


        for (int i = 0; i < numCPUs; i++){
            tasks[i] = new ParticleHistoryTask(
                    nuclearReaction, testParticleDistribution, plasma, stoppingPower, totalPerThread);

            //tasks[i].setDebugMode(true);
            tasks[i].run();

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
