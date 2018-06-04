package MonteCarloParticleTracer;

import org.apache.commons.math3.analysis.interpolation.BicubicInterpolatingFunction;
import org.apache.commons.math3.analysis.interpolation.BicubicInterpolator;

/**
 * Created by lahmann on 2017-03-01.
 */
public class StoppingPowerModel {

    private static final double MINIMUM_LI_PETRASSO_ENERGY_MeV = 0.01;

    private final double MAX_ENERGY = 20.0;
    private final double MIN_ENERGY = MINIMUM_LI_PETRASSO_ENERGY_MeV;       // MeV  TODO: Handle this better
    private final int NUM_NODES = 200;

    private BicubicInterpolatingFunction function;     // Stopping power profile defined in f(E,r) space (MeV / cm)
    private double[] radiusNodes;
    private double[] energyNodes;


    public StoppingPowerModel(ParticleType particleType, Plasma plasma){

        energyNodes = new double[NUM_NODES +1];
        radiusNodes = plasma.getNormalizedRadiusNodes();


        energyNodes[0] = 0.0;
        double[] temp = DoubleArray.linspace(MIN_ENERGY, MAX_ENERGY, NUM_NODES).getValues();
        for (int i = 1; i < energyNodes.length; i++){
            energyNodes[i] = temp[i-1];
        }


        double[][] dEdx = new double[energyNodes.length][radiusNodes.length];
        for (int i = 1; i < energyNodes.length; i++){
            dEdx[i] = plasma.getStoppingPower(particleType, energyNodes[i]);
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
        if (radius < radiusNodes[0])                        radius = radiusNodes[0];
        if (radius > radiusNodes[radiusNodes.length - 1])   radius = radiusNodes[radiusNodes.length - 1];

        if (energy < energyNodes[0])                        energy = energyNodes[0];
        if (energy > energyNodes[energyNodes.length - 1])   energy = energyNodes[energyNodes.length - 1];

        return function.value(energy, radius);
    }
}
