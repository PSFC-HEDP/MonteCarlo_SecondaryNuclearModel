import cStopPow.DoubleVector;
import cStopPow.StopPow_LP;
import org.apache.commons.math3.analysis.interpolation.BicubicInterpolatingFunction;
import org.apache.commons.math3.analysis.interpolation.BicubicInterpolator;
import org.apache.commons.math3.util.FastMath;

/**
 * Created by lahmann on 2017-03-01.
 */
public class StoppingPowerModel {

    final double MINIMUM_STOPPING_POWER_ENERGY = Utils.MINIMUM_LI_PETRASSO_ENERGY_MeV;       // MeV
    final int NUM_ENERGY_NODES = 200;

    private BicubicInterpolatingFunction function;     // Stopping power profile defined in f(E,r) space (MeV / cm)
    private double[] radiusNodes;
    private double[] energyNodes;


    public StoppingPowerModel(ParticleDistribution particleDistribution, Plasma plasma){

        final double maxEnergy = particleDistribution.getMaxEnergy();
        energyNodes = new double[NUM_ENERGY_NODES+1];
        radiusNodes = plasma.getRadiusNodes();


        energyNodes[0] = 0.0;
        double[] temp = Utils.linspace(MINIMUM_STOPPING_POWER_ENERGY, maxEnergy, NUM_ENERGY_NODES);
        for (int i = 1; i < energyNodes.length; i++){
            energyNodes[i] = temp[i-1];
        }


        double[][] dEdx = new double[energyNodes.length][radiusNodes.length];
        for (int i = 1; i < energyNodes.length; i++){
            dEdx[i] = plasma.getStoppingPower(particleDistribution.getType(), energyNodes[i]);
        }


        // Fill "zero energy point" with an extrapolation
        for (int i = 0; i < radiusNodes.length; i++){
            double E0 = energyNodes[0];
            double E1 = energyNodes[1];
            double E2 = energyNodes[2];

            dEdx[0][i] = dEdx[1][i] - (dEdx[2][i]-dEdx[1][i])*(E1-E0)/(E2-E1);
        }

        this.function = new BicubicInterpolator().interpolate(energyNodes, radiusNodes, dEdx);
    }

    public double evaluate(double energy, double radius){
        if (function.isValidPoint(energy, radius)){
            return function.value(energy, radius);
        }else{
            radius = FastMath.max(radius, radiusNodes[0]);
            radius = FastMath.min(radius, radiusNodes[radiusNodes.length-1]);

            energy = FastMath.max(energy, energyNodes[0]);
            energy = FastMath.min(energy, energyNodes[radiusNodes.length-1]);

            return function.value(energy, radius);
        }
    }
}
