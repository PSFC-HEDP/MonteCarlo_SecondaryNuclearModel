package MonteCarloParticleTracer;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.io.File;

/**
 * Particle time of flight detector
 */
public class PTOF_Detector {

    private final static Vector3D SPEC_SP_LOS =
            Utils.getVectorFromSpherical(1798, Math.toRadians(161), Math.toRadians(56));

    private final static Vector3D SPEC_A_LOS =
            Utils.getVectorFromSpherical(2222, Math.toRadians(116), Math.toRadians(316));

    private final static Vector3D SPEC_E_LOS =
            Utils.getVectorFromSpherical(2009, Math.toRadians(90), Math.toRadians(174));

    private final static Vector3D SPEC_NP_LOS =
            Utils.getVectorFromSpherical(2161, Math.toRadians(18), Math.toRadians(303));

    String name;

    private Vector3D position;
    private ParticleType particle;

    private PolynomialSplineFunction energySensitivity;
    private PolynomialSplineFunction instrumentResponse;
    private PolynomialFunction highEnergyExtrapolation;

    private double resistance;



    public static PTOF_Detector SPEC_SP01_DT(){
        PTOF_Detector spec = new PTOF_Detector("SPEC-SP01", ParticleType.neutron, SPEC_SP_LOS, 50.0);
        spec.setInstrumentResponse(new File("data/NTOF_IRFs/SPEC_SP01_DTn_IRF.dat"));
        spec.setEnergySensitivity(new File("data/NTOF_IRFs/SPEC_SP01_ENERGY_RESPONSE.dat"),
                new double[] {14.1, 1e-12 / 9.37e7});
        spec.setHighEnergyExtrapolation(spec.energySensitivity, 12.0);
        return spec;
    }

    public static PTOF_Detector SPEC_A01_DT(){
        PTOF_Detector spec = new PTOF_Detector("SPEC-A01", ParticleType.neutron, SPEC_A_LOS, 50.0);
        spec.setInstrumentResponse(new File("data/NTOF_IRFs/SPEC_A01_DTn_IRF.dat"));
        spec.setEnergySensitivity(new File("data/NTOF_IRFs/SPEC_A01_ENERGY_RESPONSE.dat"),
                new double[] {14.1, 1e-12 / 7.62e7});
        spec.setHighEnergyExtrapolation(spec.energySensitivity, 12.0);
        return spec;
    }

    public static PTOF_Detector SPEC_E01_DT(){
        PTOF_Detector spec = new PTOF_Detector("SPEC-E01", ParticleType.neutron, SPEC_E_LOS, 50.0);
        spec.setInstrumentResponse(new File("data/NTOF_IRFs/SPEC_E01_DTn_IRF.dat"));
        spec.setEnergySensitivity(new File("data/NTOF_IRFs/SPEC_E01_ENERGY_RESPONSE.dat"),
                new double[] {14.1, 1e-12 / 6.11e7});
        spec.setHighEnergyExtrapolation(spec.energySensitivity, 12.0);
        return spec;
    }

    public static PTOF_Detector SPEC_NP01_DT(){
        PTOF_Detector spec = new PTOF_Detector("SPEC-NP01", ParticleType.neutron, SPEC_NP_LOS, 50.0);
        spec.setInstrumentResponse(new File("data/NTOF_IRFs/SPEC_NP01_DTn_IRF.dat"));
        spec.setEnergySensitivity(new File("data/NTOF_IRFs/SPEC_NP01_ENERGY_RESPONSE.dat"),
                new double[] {14.1, 1e-12 / 15.6e7});
        spec.setHighEnergyExtrapolation(spec.energySensitivity, 12.0);
        return spec;
    }

    public PTOF_Detector(String name, ParticleType particle, double distance, double polarAngle, double azimuth, double resistance) {
        this.name = name;
        this.particle = particle;
        this.position = new Vector3D(distance, Math.toRadians(polarAngle), Math.toRadians(azimuth));
        this.resistance = resistance;
    }

    public PTOF_Detector(String name, ParticleType particle, Vector3D position, double resistance) {
        this.name = name;
        this.particle = particle;
        this.position = position;
        this.resistance = resistance;
    }

    public void setEnergySensitivity(File energySensitivity, double[] calibration) {
        try {
            double[][] data = Utils.parseCSV(energySensitivity);
            this.energySensitivity = new LinearInterpolator().interpolate(data[0], data[1]);

            // Apply calibration
            double factor = calibration[1] / this.energySensitivity.value(calibration[0]);
            for (int i = 0; i < data[1].length; i++){
                data[1][i] *= factor;
            }

            this.energySensitivity = new LinearInterpolator().interpolate(data[0], data[1]);
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }

    public void setInstrumentResponse(File instrumentResponse) {
        try {
            double[][] data = Utils.parseCSV(instrumentResponse);
            this.instrumentResponse = new LinearInterpolator().interpolate(data[0], data[1]);
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }

    public void setHighEnergyExtrapolation(PolynomialSplineFunction function, double startPoint) {

        SimpleRegression regression = new SimpleRegression();
        double[] x = function.getKnots();
        for (int i = 0; i < x.length; i++){
            if (x[i] > startPoint) {
                regression.addData(x[i], function.value(x[i]));
            }
        }

        highEnergyExtrapolation = new PolynomialFunction(new double[] {regression.getIntercept(), regression.getSlope()});
    }

    public String getName() {
        return name;
    }

    public Vector3D getPosition() {
        return position;
    }

    public Tally generateResponse(Tally data) throws Exception {

        // *****************************
        // Apply the energy sensitivity
        // *****************************

        double[] edges   = data.getBinEdges();
        double[] centers = data.getBinCenters();
        double[] N = data.getWeights();

        for (int i = 0; i < centers.length; i++){
            try {
                N[i] *= resistance * energySensitivity.value(centers[i]);
            }
            catch (OutOfRangeException e){
                N[i] *= resistance * highEnergyExtrapolation.value(centers[i]);
            }
        }


        // ************************************
        // Apply the time of flight conversion
        // ************************************

        for (int i = 0; i < edges.length; i++){
            if (edges[i] == 0)  {
                edges[i] = 1;
            }
            else {
                edges[i] = position.getNorm() * Math.sqrt(0.5 * Constants.MEV_PER_AMU * particle.getMass() / edges[i])
                        / Constants.SPEED_OF_LIGHT_CM_PER_SEC;
                N[i-1] /= Math.abs(edges[i] - edges[i-1]);
            }

        }


        // *****************
        // Convolve the IRF
        // *****************

        Tally tally = new Tally(edges);
        tally.weights = N;
        convolveWithIRF(tally);

        return tally;
    }

    private void convolveWithIRF(Tally data){

        double[] times  = data.getBinCenters();
        double[] raw    = data.getWeights();
        double[] signal = new double[data.getWeights().length];

        for (int i = 0; i < signal.length; i++){
            for (int j = 1; j < signal.length; j++) {

                double value0 = 0;
                try {
                    value0 = raw[j-1] * instrumentResponse.value(times[i] - times[j-1]);
                } catch (OutOfRangeException e){ }

                double value1 = 0;
                try {
                    value0 = raw[j] * instrumentResponse.value(times[i] - times[j]);
                } catch (OutOfRangeException e){ }

                signal[i] += Math.abs(times[j] - times[j-1]) * (value1 + value0);
            }
        }

        data.weights = signal;
    }
}
