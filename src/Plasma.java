import cStopPow.DoubleVector;
import cStopPow.StopPow;
import cStopPow.StopPow_LP;
import cStopPow.StopPow_Zimmerman;
import org.apache.commons.math3.analysis.function.Cos;
import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator;
import org.apache.commons.math3.analysis.integration.UnivariateIntegrator;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialsUtils;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.ArithmeticUtils;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathUtils;

import java.lang.reflect.Method;
import java.util.ArrayList;

/**
 * Created by lahmann on 2016-06-15.
 *
 * Class representing a plasma object in the model.
 * Any model should in theory only be working with one plasma for the time being.
 * This object contains the temperature and density profiles
 *
 * Currently this model assumes no species separation
 *
 */
public class Plasma {

    // Constants used for volumetric integration
    private final int NUM_PHI_NODES    = 100;
    private final int NUM_THETA_NODES  = 100;
    private final int NUM_RADIUS_NODES = 100;

    // Constants for handling floating point exceptions
    private final double OUT_BOUNDS_EPSILON = 1e-6;


    // ArrayList of Legendre Modes that describe the shape of this plasma
    private ArrayList<LegendreMode> legendreModes = new ArrayList<>();

    // ArrayList of Particle Species in this plasma
    private ArrayList<PlasmaSpecies> plasmaSpecies = new ArrayList<>();


    private PolynomialSplineFunction ionTemperature;        // Ion Temperature radial profile define between [0,1] (keV)
    private PolynomialSplineFunction electronTemperature;   // Electron Temperature radial profile defined between [0, 1] (keV)
    private PolynomialSplineFunction massDensity;           // Mass density radial profile defined between [0, 1] (g/cc)



    /**
     * Built in pre-defined plasma profiles
     */

    // Profile defined by T(r) = T and \rho(r) = \rho
    static Plasma uniformPlasma(double P0, double totalMass, double burnT){
        return uniformPlasma(P0, totalMass, burnT, 1.0);
    }

    static Plasma uniformPlasma(double P0, double totalMass, double burnT, double TeFraction){
        double[] uniformArray = Utils.linspace(1, 1, 3);

        Plasma plasma = new Plasma(uniformArray, uniformArray, uniformArray);
        plasma.addLegendreMode(0, 0, Math.sqrt(4*Math.PI)*P0);
        plasma.setTotalMass(totalMass);
        plasma.setDDnBurnAveragedIonTemperature(burnT);
        plasma.setElectronTemperatureFraction(TeFraction);
        return plasma;
    }


    // Profile defined by T(r) = T0 [1 - (r/R)^2]^{0.4} / [1 - 0.15(r/R)^2]
    // Density profile is simply \rho(r) = \rho0 (T0 / T(r))
    static Plasma BettiPlasma(double P0, double totalMass, double burnT){
        return BettiPlasma(P0, totalMass, burnT, 1.0);
    }

    static Plasma BettiPlasma(double P0, double totalMass, double burnT, double TeFraction){
        final double boundaryEpsilon  =       1e-3;     // Density profiles explode at r = R
        final  int   numProfilePoints = (int) 2e+2;

        double[] r   = Utils.linspace(0, (1-boundaryEpsilon), numProfilePoints);    // Normalized r
        double[] T   = new double[r.length];                                           // Normalized T
        double[] rho = new double[r.length];                                           // Normalized rho

        for (int i = 0; i < r.length; i++){
            T[i] = Math.pow(1 - Math.pow(r[i], 2), 0.4);
            T[i] /= (1 - 0.15*Math.pow(r[i], 2));
            rho[i] = 1.0 / T[i];
        }

        Plasma plasma = new Plasma(T, T, rho);
        plasma.addLegendreMode(0, 0, Math.sqrt(4*Math.PI)*P0);
        plasma.setTotalMass(totalMass);
        plasma.setDDnBurnAveragedIonTemperature(burnT);
        plasma.setElectronTemperatureFraction(TeFraction);
        return plasma;
    }


    // Profile defined by T(r) = T0 [1 - (r/R)^2]^{0.331}
    // Density profile is simply \rho(r) = \rho0 (T0 / T(r))
    static Plasma PravPlasma(double P0, double totalMass, double burnT){
        return PravPlasma(P0, totalMass, burnT, 1.0);
    }

    static Plasma PravPlasma(double P0, double totalMass, double burnT, double TeFraction){
        final double boundaryEpsilon  =       1e-3;     // Density profiles explode at r = R
        final  int   numProfilePoints = (int) 2e+2;

        double[] r   = Utils.linspace(0, (1-boundaryEpsilon), numProfilePoints);    // Normalized r
        double[] T   = new double[r.length];                                           // Normalized T
        double[] rho = new double[r.length];                                           // Normalized rho

        for (int i = 0; i < r.length; i++){
            T[i] = Math.pow(1 - Math.pow(r[i], 2), 0.331);
            rho[i] = 1.0 / T[i];
        }

        Plasma plasma = new Plasma(T, T, rho);
        plasma.addLegendreMode(0, 0, Math.sqrt(4*Math.PI)*P0);
        plasma.setTotalMass(totalMass);
        plasma.setDDnBurnAveragedIonTemperature(burnT);
        plasma.setElectronTemperatureFraction(TeFraction);
        return plasma;
    }


    /**
     * Constructor methods
     */

    Plasma(double[] ionTemperature, double[] electronTemperature, double[] massDensity) {
        double[] r = Utils.linspace(0, 1, ionTemperature.length);

        this.ionTemperature = new SplineInterpolator().interpolate(r, ionTemperature);
        this.electronTemperature = new SplineInterpolator().interpolate(r, electronTemperature);
        this.massDensity = new SplineInterpolator().interpolate(r, massDensity);
    }

    Plasma(double P0, double[] ionTemperature, double[] electronTemperature, double[] massDensity){
        this(ionTemperature, electronTemperature, massDensity);
        addLegendreMode(0, 0, Math.sqrt(4*Math.PI)*P0);
    }

    void addLegendreMode(int l, int m, double magnitude){
        legendreModes.add(new LegendreMode(l, m, magnitude));
    }

    void addSpecies(int Z, int A, double numberFraction){
        plasmaSpecies.add(new PlasmaSpecies(A, Z, numberFraction));
    }

    void addDeuteriumSpecies(double numberFraction){
        this.addSpecies(1, 2, numberFraction);
    }

    void addHelium3Species(double numberFraction){
        this.addSpecies(2, 3, numberFraction);
    }

    void addTritiumSpecies(double numberFraction){
        this.addSpecies(1, 3, numberFraction);
    }



    /**
     * Method that returns an array of values corresponding the the stopping power (dE/dx) of the testParticle
     * at every radius normalized location in this Plasma in units of MeV / cm
     *
     * These values are extremely sensitive to the energy of the test particle and will need to be regenerated
     * as the particle looses energy in this Plasma
     */
    double[] getStoppingPower(Particle testParticle){

        double[] r = getRadiusNodes();                  // Normalized radius [0, 1]
        double[] dEdx = new double[r.length];           // Stopping power array

        /**
         * Grab the species As and Zs since these will not very with radius
         */
        DoubleVector speciesZs = new DoubleVector();
        DoubleVector speciesAs = new DoubleVector();

        for (PlasmaSpecies species : plasmaSpecies){
            speciesZs.add(species.Z);
            speciesAs.add(species.A);
        }


        for (int i = 0; i < dEdx.length; i++){
            double rho = massDensity.value(r[i]);           // Mass density in g/cc
            double Ti  = ionTemperature.value(r[i]);        // Ion Temperature in keV
            double Te  = electronTemperature.value(r[i]);   // Electron Temperature in keV
            double n   = numberDensityFromRho(rho);

            DoubleVector speciesNs = new DoubleVector();
            DoubleVector speciesTs = new DoubleVector();
            DoubleVector speciesZbars = new DoubleVector();
            for (PlasmaSpecies species : plasmaSpecies){
                double fi = species.numberFraction;
                double ni = n * fi;

                // TODO: This chunk of code is assuming Ti is equal for all species
                double Zbar = 20*Math.sqrt(Ti);     // Using approximation for ionization from Drake
                Zbar = Math.min(Zbar, species.Z);   // Zbar cannot exceed Z

                speciesNs.add(ni);
                speciesTs.add(Ti);
                speciesZbars.add(Zbar);
            }

            StopPow_Zimmerman stopPow_zimmerman = new StopPow_Zimmerman(testParticle.getA(), testParticle.getZ(),
                    speciesAs, speciesZs, speciesTs, speciesNs, speciesZbars, Te);
            dEdx[i] = 1e4 * stopPow_zimmerman.dEdx(testParticle.getE());


            // Li-Petrasso requires we add the electrons explicitly
            //speciesZs.add(-1.0);
            //speciesAs.add(1.0 / 1836.0);
            //speciesTs.add(Te);
            //speciesNs.add(ne);

            //StopPow_LP stopPow_lp = new StopPow_LP(testParticle.A, testParticle.Z,
            //        speciesAs_withElectrons, speciesZs_withElectrons, speciesTs, speciesNs);
            //dEdx[i] = 1e4 * stopPow_lp.dEdx(testParticle.E);
        }

        return dEdx;
    }

    double[] getRadiusNodes(){
        return massDensity.getKnots();
    }



    /**
     * Methods for evaluating this Plasma's conditions at some position vector r
     */

    public double getRadiusBound(double theta, double phi){
        double radius = 0.0;
        for (LegendreMode mode : legendreModes){
            radius += mode.evaluate(theta, phi);
        }
        return radius;
    }

    public double getIonTemperature(Vector3D position){
        SphericalCoordinates coordinates = new SphericalCoordinates(position);

        // Apache is in Math notation whilst we're in Physics notation
        double theta = coordinates.getPhi();
        double phi = coordinates.getTheta();

        double maxRadius = getRadiusBound(theta, phi);
        double r = coordinates.getR() / maxRadius;

        try {
            return ionTemperature.value(r);
        }

        catch (OutOfRangeException e){
            double[] knots = ionTemperature.getKnots();
            double maxValue = knots[knots.length-1];

            if (r - maxValue < OUT_BOUNDS_EPSILON){
                return ionTemperature.value(maxValue);
            }

            System.err.printf("Warning! Mass density evaluation %.6e cm out of bounds\n", coordinates.getR() - maxRadius);
        }

        return Double.NaN;
    }

    public double getElectronTemperature(Vector3D position) {
        SphericalCoordinates coordinates = new SphericalCoordinates(position);

        // Apache is in Math notation whilst we're in Physics notation
        double theta = coordinates.getPhi();
        double phi = coordinates.getTheta();

        double maxRadius = getRadiusBound(theta, phi);
        double r = coordinates.getR() / maxRadius;

        try {
            return electronTemperature.value(r);
        }

        catch (OutOfRangeException e){
            double[] knots = electronTemperature.getKnots();
            double maxValue = knots[knots.length-1];

            if (r - maxValue < OUT_BOUNDS_EPSILON){
                return electronTemperature.value(maxValue);
            }

            System.err.printf("Warning! Mass density evaluation %.6e cm out of bounds\n", coordinates.getR() - maxRadius);
        }

        return Double.NaN;
    }

    public double getDDnReactivity(Vector3D r){
        return Reactivity.DDn_Reactivity.evaluate(getIonTemperature(r));
    }

    public double getMassDensity(Vector3D position){
        SphericalCoordinates coordinates = new SphericalCoordinates(position);

        // Apache is in Math notation whilst we're in Physics notation
        double theta = coordinates.getPhi();
        double phi = coordinates.getTheta();

        double maxRadius = getRadiusBound(theta, phi);
        double r = coordinates.getR() / maxRadius;

        try {
            return massDensity.value(r);
        }

        catch (OutOfRangeException e){
            double[] knots = massDensity.getKnots();
            double maxValue = knots[knots.length-1];

            if (r - maxValue < OUT_BOUNDS_EPSILON){
                return massDensity.value(maxValue);
            }

            System.err.printf("Warning! Mass density evaluation %.6e cm out of bounds\n", coordinates.getR() - maxRadius);
        }

        return Double.NaN;
    }

    public double getNumberDensity(Vector3D r){
        return numberDensityFromRho(getMassDensity(r));
    }

    public double getDeuteronNumberDensity(Vector3D r){

        double totalFraction = 0.0;         // Sum of all species number fractions (should be 1 nominally)
        double deuteriumFraction = 0.0;     // Number fraction of deuterium atoms

        for (PlasmaSpecies species : plasmaSpecies){
            totalFraction += species.numberFraction;

            if (species.A == 2 && species.Z == 1){
                deuteriumFraction = species.numberFraction;
            }
        }

        return getNumberDensity(r) * deuteriumFraction / totalFraction;
    }

    public boolean getIsInside(Vector3D r){
        SphericalCoordinates coordinates = new SphericalCoordinates(r);

        // Apache is in Math notation whilst we're in Physics notation
        double theta = coordinates.getPhi();
        double phi = coordinates.getTheta();

        return r.getNorm() < getRadiusBound(theta, phi);
    }


    /**
     * Methods for evaluating volume integrated quantities of this Plasma
     */

    public double getTotalVolume(){
        return volumeIntegral(new String[] {});
    }

    public double getTotalMass(){
        return volumeIntegral("getMassDensity");
    }

    public double getVolumeAverageIonTemperature(){
        return volumeIntegral("getIonTemperature") / getTotalVolume();
    }

    public double getVolumeAverageElectronTemperature(){
        return volumeIntegral("getElectronTemperature") / getTotalVolume();
    }

    // Get average Tion weighted by rho^2 <sigmaV>
    public double getDDnBurnAveragedIonTemperature(){
        double numerator   = volumeIntegral("getMassDensity", "getMassDensity", "getDDnReactivity", "getIonTemperature");
        double denominator = volumeIntegral("getMassDensity", "getMassDensity", "getDDnReactivity");
        return  numerator / denominator;
    }


    /**
     * Convenience method for generating spatial burn distributions using this Plasma's profiles
     * Distributions are defined in normalized radius [0, 1]
     */
    Distribution getSpatialDDnBurnDistribution(){
        return getSpatialBurnDistribution(Reactivity.DDn_Reactivity);
    }


    Distribution getSpatialDDpBurnDistribution(){
        return getSpatialBurnDistribution(Reactivity.DDp_Reactivity);
    }

    Distribution getSpatialBurnDistribution(Reactivity reactivity){
        double[] r = getRadiusNodes();                  // Normalized r
        double[] Y = new double[r.length];              // Yield probability

        for (int i = 0; i < r.length; i++){
            double rho = massDensity.value(r[i]);
            double T   = ionTemperature.value(r[i]);
            double sigmaV = reactivity.evaluate(T);

            Y[i] = r[i] * r[i] * rho * rho * sigmaV;        // Y(r) \proto n^2 <\sigma v>. r^2 factor because spherical
        }

        return new Distribution(r, Y);
    }



    /**
     * Methods for setting volume integrated quantities of this Plasma
     */

    public void setTotalMass(double totalMass){
        double multiplicativeFactor = totalMass / getTotalMass();
        PolynomialFunction constantPoly = new PolynomialFunction(new double[] {multiplicativeFactor});

        PolynomialFunction[] polynomialFunctions = massDensity.getPolynomials();
        for (int i = 0; i < polynomialFunctions.length; i++){
            polynomialFunctions[i] = polynomialFunctions[i].multiply(constantPoly);
        }

        massDensity = new PolynomialSplineFunction(massDensity.getKnots(), polynomialFunctions);
    }

    public void setDDnBurnAveragedIonTemperature(double burnT){
        double multiplicativeFactor = burnT / getDDnBurnAveragedIonTemperature();
        PolynomialFunction constantPoly = new PolynomialFunction(new double[] {multiplicativeFactor});

        PolynomialFunction[] polynomialFunctions = ionTemperature.getPolynomials();
        for (int i = 0; i < polynomialFunctions.length; i++){
            polynomialFunctions[i] = polynomialFunctions[i].multiply(constantPoly);
        }

        ionTemperature = new PolynomialSplineFunction(ionTemperature.getKnots(), polynomialFunctions);
    }

    public void setElectronTemperatureFraction(double electronTemperatureFraction){
        double multiplicativeFactor = electronTemperatureFraction
                * getVolumeAverageIonTemperature() / getVolumeAverageElectronTemperature();
        PolynomialFunction constantPoly = new PolynomialFunction(new double[] {multiplicativeFactor});

        PolynomialFunction[] polynomialFunctions = electronTemperature.getPolynomials();
        for (int i = 0; i < polynomialFunctions.length; i++){
            polynomialFunctions[i] = polynomialFunctions[i].multiply(constantPoly);
        }

        electronTemperature = new PolynomialSplineFunction(electronTemperature.getKnots(), polynomialFunctions);
    }



    /**
     * Class for Legendre Mode information
     */
    class LegendreMode{
        private int el, m;
        private double norm;
        private double magnitude;       // um

        LegendreMode(int el, int m, double magnitude) {
            this.el = el;
            this.m = m;
            this.magnitude = magnitude;

            this.norm = CombinatoricsUtils.factorial(el - m);
            norm /= CombinatoricsUtils.factorial(el + m);
            norm *= (2*el + 1);
            norm /= (4 * Math.PI);
            norm = Math.sqrt(norm);
        }

        double evaluate(double theta, double phi){
            double value = this.magnitude;
            value *= this.norm;
            value *= FastMath.cos(m * phi);
            value *= associatedLegendrePolyValue(theta);
            return value;
        }

        private double associatedLegendrePolyValue(double theta){
            double x = Math.cos(theta);     // Evaluation point
            double p0 = 1.0;
            int el = 0, m = 0;

            // Recurse until we find P(x, el=M, m=M)
            while ( m < this.m){
                p0 *= -(2*el + 1) * Math.sin(theta);
                m++;
                el++;
            }

            // If that's the one we need, return it
            if (el == this.el) return p0;

            // Else calculate P(x, el=M+1, m=M)
            double p1 = p0 * x * (2*el + 1);
            el++;

            // Recurse until we fine P(x, el=L, m=M)
            while ( el < this.el){
                double temp = p1;

                p1 = (2*el + 1)*x*p1 - (el + m)*p0;
                p1 /= (el - m + 1);

                p0 = temp;
                el++;
            }

            return p1;
        }
    }


    /**
     * Glorified struct to hold Plasma Species information
     */
    class PlasmaSpecies {

        int A, Z;
        double numberFraction;

        PlasmaSpecies(int a, int z, double numberFraction) {
            A = a;
            Z = z;
            this.numberFraction = numberFraction;
        }
    }


    /**
     * Private convenience methods
     */

    private double getAbar(){
        double Abar = 0.0;
        double totalFraction = 0.0;

        for (PlasmaSpecies species : plasmaSpecies){
            Abar += species.A * species.numberFraction;
            totalFraction += species.numberFraction;
        }

        return Abar / totalFraction;
    }

    private double numberDensityFromRho(double rho){
        return rho / (getAbar() * 1000 * Utils.protonMass_kg);
    }

    private double volumeIntegral(String ... methodNames){
        Method[] methods = new Method[methodNames.length];

        try {

            for (int i = 0; i < methods.length; i++){
                methods[i] = this.getClass().getMethod(methodNames[i], Vector3D.class);
            }

        }catch (NoSuchMethodException e){
            e.printStackTrace();
            System.exit(-1);
        }

        return volumeIntegral(methods);
    }

    private double volumeIntegral(Method ... weightFunctions){
        try {
            final double MAX_PHI = 2 * Math.PI;
            final double MAX_THETA = 1 * Math.PI;

            double total = 0.0;
            for (int k = 0; k < NUM_PHI_NODES; k++) {

                double phi_1 = (k + 1) * MAX_PHI / (NUM_PHI_NODES);
                double phi_0 = (k + 0) * MAX_PHI / (NUM_PHI_NODES);
                double dPhi = (phi_1 - phi_0);

                for (int j = 0; j < NUM_THETA_NODES; j++) {

                    double theta_1 = (j + 1) * MAX_THETA / (NUM_THETA_NODES);
                    double theta_0 = (j + 0) * MAX_THETA / (NUM_THETA_NODES);
                    double dTheta = (theta_1 - theta_0);

                    for (int i = 0; i < NUM_RADIUS_NODES; i++) {

                        double total_1 = 0.0;       // Half of the total corresponding to theta_1
                        double total_0 = 0.0;       // Half of the total corresponding to theta_0

                        double tempHalf_1 = 0.0;    // Reusable temp product corresponding to r_1
                        double tempHalf_0 = 0.0;    // Reusable temp product corresponding to r_0

                        // All theta_1 and phi_1 quantities
                        double r_111 = (i + 1) * getRadiusBound(theta_1, phi_1) / NUM_RADIUS_NODES;
                        double r_011 = (i + 0) * getRadiusBound(theta_1, phi_1) / NUM_RADIUS_NODES;
                        Vector3D vec_111 = new SphericalCoordinates(r_111, theta_1, phi_1).getCartesian();
                        Vector3D vec_011 = new SphericalCoordinates(r_011, theta_1, phi_1).getCartesian();
                        double dR_11 = (r_111 - r_011);


                        // Add them to total_1
                        tempHalf_1 = r_111 * r_111;
                        for (Method weight : weightFunctions){
                            tempHalf_1 *= (double) weight.invoke(this, vec_111);
                        }

                        tempHalf_0 = r_011 * r_011;
                        for (Method weight : weightFunctions){
                            tempHalf_0 *= (double) weight.invoke(this, vec_011);
                        }

                        total_1 += dR_11 * (tempHalf_1 + tempHalf_0);


                        // All theta_1 and phi_0 quantities
                        double r_110 = (i + 1) * getRadiusBound(theta_1, phi_0) / NUM_RADIUS_NODES;
                        double r_010 = (i + 0) * getRadiusBound(theta_1, phi_0) / NUM_RADIUS_NODES;
                        Vector3D vec_110 = new SphericalCoordinates(r_110, theta_1, phi_1).getCartesian();
                        Vector3D vec_010 = new SphericalCoordinates(r_010, theta_1, phi_1).getCartesian();
                        double dR_10 = (r_110 - r_010);


                        // Add them to total_1
                        tempHalf_1 = r_110 * r_110;
                        for (Method weight : weightFunctions){
                            tempHalf_1 *= (double) weight.invoke(this, vec_110);
                        }

                        tempHalf_0 = r_010 * r_010;
                        for (Method weight : weightFunctions){
                            tempHalf_0 *= (double) weight.invoke(this, vec_010);
                        }

                        total_1 += dR_10 * (tempHalf_1 + tempHalf_0);


                        // All theta_0 and phi_1 quantities
                        double r_101 = (i + 1) * getRadiusBound(theta_0, phi_1) / NUM_RADIUS_NODES;
                        double r_001 = (i + 0) * getRadiusBound(theta_0, phi_1) / NUM_RADIUS_NODES;
                        Vector3D vec_101 = new SphericalCoordinates(r_101, theta_1, phi_1).getCartesian();
                        Vector3D vec_001 = new SphericalCoordinates(r_001, theta_1, phi_1).getCartesian();
                        double dR_01 = (r_101 - r_001);


                        // Add them to total_0
                        tempHalf_1 = r_101 * r_101;
                        for (Method weight : weightFunctions){
                            tempHalf_1 *= (double) weight.invoke(this, vec_101);
                        }

                        tempHalf_0 = r_001 * r_001;
                        for (Method weight : weightFunctions){
                            tempHalf_0 *= (double) weight.invoke(this, vec_001);
                        }

                        total_0 += dR_01 * (tempHalf_1 + tempHalf_0);


                        // All theta_0 and phi_0 quantities
                        double r_100 = (i + 1) * getRadiusBound(theta_0, phi_0) / NUM_RADIUS_NODES;
                        double r_000 = (i + 0) * getRadiusBound(theta_0, phi_0) / NUM_RADIUS_NODES;
                        Vector3D vec_100 = new SphericalCoordinates(r_100, theta_1, phi_1).getCartesian();
                        Vector3D vec_000 = new SphericalCoordinates(r_000, theta_1, phi_1).getCartesian();
                        double dR_00 = (r_100 - r_000);


                        // Add them to total_0
                        tempHalf_1 = r_100 * r_100;
                        for (Method weight : weightFunctions){
                            tempHalf_1 *= (double) weight.invoke(this, vec_100);
                        }

                        tempHalf_0 = r_000 * r_000;
                        for (Method weight : weightFunctions){
                            tempHalf_0 *= (double) weight.invoke(this, vec_000);
                        }

                        total_0 += dR_00 * (tempHalf_1 + tempHalf_0);


                        // Finally add them together with appropriate weights
                        total += (dPhi * dTheta / 8) * (FastMath.sin(theta_1) * total_1 + FastMath.sin(theta_0) * total_0);
                    }
                }
            }

            return total;
        }catch (Exception e){
            e.printStackTrace();
            System.exit(-1);
        }

        return Double.NaN;
    }

}