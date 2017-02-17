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

    // Constants used when generating the stopping power function
    final double MINIMUM_STOPPING_POWER_ENERGY = 0.1;       // MeV
    final int NUM_ENERGY_NODES = 50;

    private ParticleDistribution testParticleDistribution;
    private BicubicInterpolatingFunction stoppingPower;     // Stopping power profile defined in f(E,r) space (MeV / cm)
    private CrossSection crossSection;
    private Plasma plasma;



    public Model(ParticleDistribution testParticleDistribution, CrossSection crossSection, Plasma plasma) {
        this.testParticleDistribution = testParticleDistribution;
        this.crossSection = crossSection;
        this.plasma = plasma;
        this.stoppingPower = calculateStoppingPowerFunction();
    }

    public double getYieldRatio(int totalParticles){

        //int numCPUs = Runtime.getRuntime().availableProcessors();
        int numCPUs = 1;
        int totalPerThread = (int) Math.ceil((double) totalParticles / numCPUs);

        Thread[] threads = new Thread[numCPUs];
        ParticleHistoryTask[] tasks = new ParticleHistoryTask[numCPUs];


        for (int i = 0; i < numCPUs; i++){
            tasks[i] = new ParticleHistoryTask(
                    crossSection, testParticleDistribution, plasma, stoppingPower, totalPerThread);
            //tasks[i].setDebugMode(true);
            threads[i] = new Thread(tasks[i]);
            threads[i].start();
        }

        double yieldRatio = 0.0;
        for (int i = 0; i < numCPUs; i ++){
            try {
                threads[i].join();
                yieldRatio += (double) tasks[i].get();

            }catch (Exception e){
                e.printStackTrace();
            }
        }

        return yieldRatio/numCPUs;
    }

    private BicubicInterpolatingFunction calculateStoppingPowerFunction(){

        final double maxEnergy = testParticleDistribution.getMaxEnergy();
        double[] energyNodes = Utils.linspace(MINIMUM_STOPPING_POWER_ENERGY, maxEnergy, NUM_ENERGY_NODES);
        double[] radiusNodes = plasma.getRadiusNodes();

        double[][] dEdx = new double[energyNodes.length][radiusNodes.length];
        for (int i = 0; i < energyNodes.length; i++){
            Particle p = new Particle(energyNodes[i], testParticleDistribution.getZ(), testParticleDistribution.getA());
            dEdx[i] = plasma.getStoppingPower(p);
        }

        return new BicubicInterpolator().interpolate(energyNodes, radiusNodes, dEdx);
    }

}
