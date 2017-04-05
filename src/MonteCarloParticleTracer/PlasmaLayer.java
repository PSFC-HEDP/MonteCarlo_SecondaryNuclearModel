package MonteCarloParticleTracer;

import cStopPow.DoubleVector;
import cStopPow.StopPow_LP;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;

import java.lang.reflect.Method;
import java.util.ArrayList;

/**
 * This class represents a single layer of plasma defined by radial profiles of mass density, ion temperature, and
 * electron temperature.
 * The plasma is bounded by two arbitrary boundaries defined by Legendre modes.
 * The plasma profiles are assumed to be separable in r, theta, and phi such that the boundaries and radial profiles
 * can be stored as separate Objects
 *
 * PlasmaLayers are "stacked" in the MonteCarloParticleTracer.Plasma Object to represent distinct layers
 *
 * Currently this model assumes no species separation (number fractions are scalars as opposed to profiles)
 *
 */
public class PlasmaLayer {

    // Constants used for volumetric integration
    // Benchmark note (50/50/50 gets sphere volume correct with 0.013%)
    private static final int NUM_PHI_NODES    = 50;
    private static final int NUM_THETA_NODES  = 50;
    private static final int NUM_RADIUS_NODES = 50;

    // Constants for handling floating point exceptions
    private static final double OUT_BOUNDS_EPSILON = 1e-6;

    // Constants for predefined profiles
    private static final int NUM_PROFILE_POINTS = 200;


    // ArrayList of Legendre Modes that describe the inner boundary of this PlasmaLayer
    private ArrayList<LegendreMode> innerBoundaryLegendreModes = new ArrayList<>();

    // ArrayList of Legendre Modes that describe the outer boundary of this PlasmaLayer
    private ArrayList<LegendreMode> outerBoundaryLegendreModes = new ArrayList<>();

    // ArrayList of Particle Species in this PlasmaLayer
    private ArrayList<PlasmaSpecies> plasmaSpecies = new ArrayList<>();

    // Ion Temperature radial profile define between [rMin, rMax] (keV)
    // Internally, nodes are mapped to the normalized units of [0, 1]
    private PolynomialSplineFunction ionTemperature;

    // Electron Temperature radial profile defined between [rMin, rMax] (keV)
    // Internally, nodes are mapped to the normalized units of [0, 1]
    private PolynomialSplineFunction electronTemperature;

    // Mass density radial profile defined between [rMin, rMax] (keV)
    // Internally, nodes are mapped to the normalized units of [0, 1]
    private PolynomialSplineFunction massDensity;



    /**
     * Built in pre-defined plasma profiles
     */

    // Profile defined by T(r) = T and \rho(r) = \rho
    public static PlasmaLayer uniformPlasma(double P0, double totalMass, double burnT){
        return uniformPlasma(P0, totalMass, burnT, 1.0);
    }

    public static PlasmaLayer uniformPlasma(double P0, double totalMass, double burnT, double TeFraction){
        double[] uniformArray = Utils.linspace(1, 1, NUM_PROFILE_POINTS);

        PlasmaLayer plasma = new PlasmaLayer(uniformArray, uniformArray, uniformArray);
        plasma.setOuterP0(P0);
        plasma.setTotalMass(totalMass);
        plasma.setDDnBurnAveragedIonTemperature(burnT);
        plasma.setElectronTemperatureFraction(TeFraction);
        return plasma;
    }


    // Profile defined by T(r) = T0 [1 - (r/R)^2]^{0.4} / [1 - 0.15(r/R)^2]
    // Density profile is simply \rho(r) = \rho0 (T0 / T(r))
    public static PlasmaLayer BettiPlasma(double P0, double totalMass, double burnT){
        return BettiPlasma(P0, totalMass, burnT, 1.0);
    }

    public static PlasmaLayer BettiPlasma(double P0, double totalMass, double burnT, double TeFraction){
        final double boundaryEpsilon  =       1e-3;     // Density profiles explode at r = R

        double[] r   = Utils.linspace(0, (1-boundaryEpsilon), NUM_PROFILE_POINTS);  // Normalized r
        double[] T   = new double[r.length];                                           // Normalized T
        double[] rho = new double[r.length];                                           // Normalized rho

        for (int i = 0; i < r.length; i++){
            T[i] = Math.pow(1 - Math.pow(r[i], 2), 0.4);
            T[i] /= (1 - 0.15*Math.pow(r[i], 2));
            rho[i] = 1.0 / T[i];
        }

        PlasmaLayer plasma = new PlasmaLayer(T, T, rho);
        plasma.setOuterP0(P0);
        plasma.setTotalMass(totalMass);
        plasma.setDDnBurnAveragedIonTemperature(burnT);
        plasma.setElectronTemperatureFraction(TeFraction);
        return plasma;
    }


    // Profile defined by T(r) = T0 [1 - (r/R)^2]^{0.331}
    // Density profile is simply \rho(r) = \rho0 (T0 / T(r))
    public static PlasmaLayer PravPlasma(double P0, double totalMass, double burnT){
        return PravPlasma(P0, totalMass, burnT, 1.0);
    }

    public static PlasmaLayer PravPlasma(double P0, double totalMass, double burnT, double TeFraction){
        final double boundaryEpsilon  =       1e-3;     // Density profiles explode at r = R

        double[] r   = Utils.linspace(0, (1-boundaryEpsilon), NUM_PROFILE_POINTS);  // Normalized r
        double[] T   = new double[r.length];                                           // Normalized T
        double[] rho = new double[r.length];                                           // Normalized rho

        for (int i = 0; i < r.length; i++){
            T[i] = Math.pow(1 - Math.pow(r[i], 2), 0.331);
            rho[i] = 1.0 / T[i];
        }

        PlasmaLayer plasma = new PlasmaLayer(T, T, rho);
        plasma.setOuterP0(P0);
        plasma.setTotalMass(totalMass);
        plasma.setDDnBurnAveragedIonTemperature(burnT);
        plasma.setElectronTemperatureFraction(TeFraction);
        return plasma;
    }


    /**
     * Generic constructor methods
     */

    public PlasmaLayer(double[] ionTemperature, double[] electronTemperature, double[] massDensity) {
        double[] r = Utils.linspace(0, 1, ionTemperature.length);

        this.ionTemperature = new SplineInterpolator().interpolate(r, ionTemperature);
        this.electronTemperature = new SplineInterpolator().interpolate(r, electronTemperature);
        this.massDensity = new SplineInterpolator().interpolate(r, massDensity);
    }

    public PlasmaLayer(double P0, double[] ionTemperature, double[] electronTemperature, double[] massDensity){
        this(ionTemperature, electronTemperature, massDensity);
        setOuterP0(P0);
    }

    public void addInnerLegendreMode(int l, int m, double magnitude){
        innerBoundaryLegendreModes.add(new LegendreMode(l, m, magnitude));
    }

    public void addOuterLegendreMode(int l, int m, double magnitude){
        outerBoundaryLegendreModes.add(new LegendreMode(l, m, magnitude));
    }

    public void addSpecies(ParticleType type, double numberFraction){
        plasmaSpecies.add(new PlasmaSpecies(type, numberFraction));
    }

    public void addDeuteriumSpecies(double numberFraction){
        this.addSpecies(ParticleType.deuteron, numberFraction);
    }

    public void addHelium3Species(double numberFraction){
        this.addSpecies(ParticleType.helium3, numberFraction);
    }

    public void addTritiumSpecies(double numberFraction){
        this.addSpecies(ParticleType.triton, numberFraction);
    }

    public boolean containsSpecies(ParticleType type){

        for (PlasmaSpecies species : plasmaSpecies){
            if (species.getType().equals(type))  return true;
        }

        return false;
    }

    /**
     * Method that returns an array of values corresponding the the stopping power (dE/dx) of the testParticle
     * at every radius normalized location in this PlasmaLayer in units of MeV / cm
     *
     * These values are extremely sensitive to the energy of the test particle and will need to be regenerated
     * as the particle looses energy in this PlasmaLayer
     */
    double[] getStoppingPower(ParticleType testParticle, double energy){

        double[] r = getRadiusNodes();                  // Normalized radius [0, 1]
        double[] dEdx = new double[r.length];           // Stopping power array


        for (int i = 0; i < dEdx.length; i++){

            // Grab everything that we'll need at this location

            double rho = massDensity.value(r[i]);           // Mass density in g/cc
            double Ti  = ionTemperature.value(r[i]);        // Ion Temperature in keV
            double Te  = electronTemperature.value(r[i]);   // Electron Temperature in keV
            double n   = numberDensityFromRho(rho);         // Number density

            // Initialize the vectors for the stopping power model
            DoubleVector speciesZs = new DoubleVector();
            DoubleVector speciesZbars = new DoubleVector();
            DoubleVector speciesAs = new DoubleVector();
            DoubleVector speciesNs = new DoubleVector();
            DoubleVector speciesTs = new DoubleVector();

            // Loop through all the plasma species
            for (PlasmaSpecies species : plasmaSpecies){
                double fi = species.getNumberFraction();    // Number fraction
                double ni = n * fi;                         // Ion Density

                double Zbar = 20*Math.sqrt(Ti);             // Using approximation for ionization from Drake
                Zbar = Math.min(Zbar, species.getZ());      // Zbar cannot exceed Z

                // Add values to the vectors
                speciesZs.add(species.getZ());
                speciesZbars.add(Zbar);
                speciesAs.add(species.getMass());
                speciesNs.add(ni);
                speciesTs.add(Ti);
            }

            //StopPow_Zimmerman stopPow_zimmerman = new StopPow_Zimmerman(testParticle.getA(), testParticle.getZ(),
            //        speciesAs, speciesZs, speciesTs, speciesNs, speciesZbars, Te);
            //dEdx[i] = 1e4 * stopPow_zimmerman.dEdx_MeV_um(testParticle.getEnergy());

            StopPow_LP stopPow_lp = new StopPow_LP(testParticle.getMass(), testParticle.getZ(),
                    speciesAs, speciesZs, speciesTs, speciesNs, Te);
            dEdx[i] = 1e4 * stopPow_lp.dEdx_MeV_um(energy);
        }

        return dEdx;
    }

    public double[] getRadiusNodes(){
        return massDensity.getKnots();
    }

    public ArrayList<LegendreMode> getInnerBoundaryLegendreModes() {
        return innerBoundaryLegendreModes;
    }

    public ArrayList<LegendreMode> getOuterBoundaryLegendreModes() {
        return outerBoundaryLegendreModes;
    }

    /**
     * Methods for evaluating this PlasmaLayer's conditions at some position vector r
     */

    public double getInnerRadiusBound(double theta, double phi){
        double radius = 0.0;
        for (LegendreMode mode : innerBoundaryLegendreModes){
            radius += mode.evaluate(theta, phi);
        }
        return radius;
    }

    public double getOuterRadiusBound(double theta, double phi){
        double radius = 0.0;
        for (LegendreMode mode : outerBoundaryLegendreModes){
            radius += mode.evaluate(theta, phi);
        }
        return radius;
    }

    public double getIonTemperature(Vector3D position){
        double[] coordinates = Utils.getSphericalFromVector(position);

        double R     = coordinates[0];
        double theta = coordinates[1];
        double phi   = coordinates[2];

        double minRadius  = getInnerRadiusBound(theta, phi);
        double maxRadius  = getOuterRadiusBound(theta, phi);
        double normRadius = (R - minRadius) / (maxRadius - minRadius);

        try {
            return ionTemperature.value(normRadius);
        }

        catch (OutOfRangeException e){
            double[] knots = ionTemperature.getKnots();
            double maxValue = knots[knots.length-1];

            if (normRadius - maxValue < OUT_BOUNDS_EPSILON){
                return ionTemperature.value(maxValue);
            }

            System.err.printf("Warning! Ion temperature evaluation %.6e cm out of bounds\n", R - maxRadius);
        }

        return Double.NaN;
    }

    public double getElectronTemperature(Vector3D position) {
        double[] coordinates = Utils.getSphericalFromVector(position);

        double R     = coordinates[0];
        double theta = coordinates[1];
        double phi   = coordinates[2];

        double minRadius  = getInnerRadiusBound(theta, phi);
        double maxRadius  = getOuterRadiusBound(theta, phi);
        double normRadius = (R - minRadius) / (maxRadius - minRadius);

        try {
            return electronTemperature.value(normRadius);
        }

        catch (OutOfRangeException e){
            double[] knots = electronTemperature.getKnots();
            double maxValue = knots[knots.length-1];

            if (normRadius - maxValue < OUT_BOUNDS_EPSILON){
                return electronTemperature.value(maxValue);
            }

            System.err.printf("Warning! Electron temperature evaluation %.6e cm out of bounds\n", R - maxRadius);
        }

        return Double.NaN;
    }

    public double getReactivity(Vector3D r, Reactivity reactivity){
        return reactivity.evaluate(getIonTemperature(r));
    }

    public double getDDnReactivity(Vector3D r){
        return getReactivity(r, Reactivity.DDn_Reactivity);
    }

    public double getD3HepReactivity(Vector3D r){
        return getReactivity(r, Reactivity.D3Hep_Reactivity);
    }

    public double getMassDensity(Vector3D position){
        double[] coordinates = Utils.getSphericalFromVector(position);

        double R     = coordinates[0];
        double theta = coordinates[1];
        double phi   = coordinates[2];

        double minRadius  = getInnerRadiusBound(theta, phi);
        double maxRadius  = getOuterRadiusBound(theta, phi);
        double normRadius = (R - minRadius) / (maxRadius - minRadius);

        try {
            return massDensity.value(normRadius);
        }

        catch (OutOfRangeException e){
            double[] knots = massDensity.getKnots();
            double maxValue = knots[knots.length-1];

            if (normRadius - maxValue < OUT_BOUNDS_EPSILON){
                return massDensity.value(maxValue);
            }

            System.err.printf("Warning! Mass density evaluation %.6e cm out of bounds\n", R - maxRadius);
        }

        return Double.NaN;
    }

    public double getNumberDensity(Vector3D r){
        return numberDensityFromRho(getMassDensity(r));
    }

    public double getSpeciesNumberDensity(Vector3D r, ParticleType type){

        double totalFraction = 0.0;         // Sum of all species number fractions (should be 1 nominally)
        double speciesFraction = 0.0;       // Number fraction of the specific species

        for (PlasmaSpecies species : plasmaSpecies){
            totalFraction += species.getNumberFraction();

            if (species.getType().equals(type)){
                speciesFraction = species.getNumberFraction();
            }
        }

        return getNumberDensity(r) * speciesFraction / totalFraction;
    }

    public boolean getIsInside(Vector3D r){
        double[] coordinates = Utils.getSphericalFromVector(r);
        double theta = coordinates[1];
        double phi   = coordinates[2];

        if (r.getNorm() > getInnerRadiusBound(theta, phi)){
            if (r.getNorm() < getOuterRadiusBound(theta, phi)){
                return true;
            }
        }

        return false;
    }

    public boolean getIsCloserToOuterBoundary(Vector3D r){
        double[] coordinates = Utils.getSphericalFromVector(r);
        double theta = coordinates[1];
        double phi   = coordinates[2];

        double distanceToOuter = getOuterRadiusBound(theta, phi) - r.getNorm();
        double distanceToInner = r.getNorm() - getInnerRadiusBound(theta, phi);

        return distanceToOuter < distanceToInner;
    }




    /**
     * Methods for evaluating volume integrated quantities of this PlasmaLayer
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
    public double getBurnAveragedIonTemperature(String reactivityMethodName){
        double numerator   = volumeIntegral("getMassDensity", "getMassDensity", reactivityMethodName, "getIonTemperature");
        double denominator = volumeIntegral("getMassDensity", "getMassDensity", reactivityMethodName);
        return  numerator / denominator;

    }

    public double getDDnBurnAveragedIonTemperature(){
        return getBurnAveragedIonTemperature("getDDnReactivity");
    }

    public double getD3HepBurnAveragedIonTemperature(){
        return getBurnAveragedIonTemperature("getD3HepReactivity");
    }


    /**
     * Convenience method for generating spatial burn distributions using this PlasmaLayer's profiles
     * Distributions are defined in normalized radius [0, 1]
     */
    public Distribution getSpatialDDnBurnDistribution(){
        return getSpatialBurnDistribution(Reactivity.DDn_Reactivity);
    }

    public Distribution getSpatialDDpBurnDistribution(){
        return getSpatialBurnDistribution(Reactivity.DDp_Reactivity);
    }

    public Distribution getSpatialD3HepBurnDistribution(){
        return getSpatialBurnDistribution(Reactivity.D3Hep_Reactivity);
    }

    public Distribution getSpatialBurnDistribution(Reactivity reactivity){
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
     * Various setter methods
     */

    public void setOuterP0(double P0){

        // The P0 accepted by this function is the P0 AFTER normalization.
        // Legendre modes are defined by magnitudes BEFORE normalization so we'll need to convert
        // This is mostly to stay consistent with literature
        double unnormalizedP0 = Math.sqrt(4*Math.PI) * P0;

        for (LegendreMode mode : outerBoundaryLegendreModes){
            if (mode.getEl() == 0 && mode.getM() == 0){
                mode.setMagnitude(unnormalizedP0);
                return;
            }
        }
        addOuterLegendreMode(0, 0, unnormalizedP0);
    }

    public void setTotalMass(double totalMass) {
        double multiplicativeFactor = totalMass / getTotalMass();
        multiplyMassDensityByScalar(multiplicativeFactor);
    }

    public void setDDnBurnAveragedIonTemperature(double burnT){
        double multiplicativeFactor = burnT / getDDnBurnAveragedIonTemperature();
        multiplyIonTemperatureByScalar(multiplicativeFactor);
    }

    public void setD3HepBurnAveragedIonTemperature(double burnT){
        double multiplicativeFactor = burnT / getD3HepBurnAveragedIonTemperature();
        multiplyIonTemperatureByScalar(multiplicativeFactor);
    }

    public void setElectronTemperatureFraction(double electronTemperatureFraction){
        double multiplicativeFactor = electronTemperatureFraction
                * getVolumeAverageIonTemperature() / getVolumeAverageElectronTemperature();
        multiplyElectronTemperatureByScalar(multiplicativeFactor);
    }

    public void setCenterMassDensity(double massDensity){
        double multiplicativeFactor = massDensity / getMassDensity(new Vector3D(0,0,0));
        multiplyMassDensityByScalar(multiplicativeFactor);
    }

    public void setCenterIonTemperature(double ionTemperature){
        double multiplicativeFactor = ionTemperature / getIonTemperature(new Vector3D(0,0,0));
        multiplyIonTemperatureByScalar(multiplicativeFactor);
    }

    public void setCenterElectronTemperature(double electronTemperature){
        double multiplicativeFactor = electronTemperature / getElectronTemperature(new Vector3D(0,0,0));
        multiplyElectronTemperatureByScalar(multiplicativeFactor);
    }

    public void setInnerBoundaryLegendreModes(ArrayList<LegendreMode> innerBoundaryLegendreModes) {
        this.innerBoundaryLegendreModes = innerBoundaryLegendreModes;
    }

    public void setOuterBoundaryLegendreModes(ArrayList<LegendreMode> outerBoundaryLegendreModes) {
        this.outerBoundaryLegendreModes = outerBoundaryLegendreModes;
    }

    /**
     * Methods for multiplying the profiles by a scalar
     * Useful for avoiding the computationally expensive volume integrals associated with making new
     * PlasmaLayer objects
     */

    void multiplyMassDensityByScalar(double scalar){
        PolynomialFunction constantPoly = new PolynomialFunction(new double[] {scalar});

        PolynomialFunction[] polynomialFunctions = massDensity.getPolynomials();
        for (int i = 0; i < polynomialFunctions.length; i++){
            polynomialFunctions[i] = polynomialFunctions[i].multiply(constantPoly);
        }

        massDensity = new PolynomialSplineFunction(massDensity.getKnots(), polynomialFunctions);
    }

    void multiplyIonTemperatureByScalar(double scalar){
        PolynomialFunction constantPoly = new PolynomialFunction(new double[] {scalar});

        PolynomialFunction[] polynomialFunctions = ionTemperature.getPolynomials();
        for (int i = 0; i < polynomialFunctions.length; i++){
            polynomialFunctions[i] = polynomialFunctions[i].multiply(constantPoly);
        }

        ionTemperature = new PolynomialSplineFunction(ionTemperature.getKnots(), polynomialFunctions);
    }

    void multiplyElectronTemperatureByScalar(double scalar){
        PolynomialFunction constantPoly = new PolynomialFunction(new double[] {scalar});

        PolynomialFunction[] polynomialFunctions = electronTemperature.getPolynomials();
        for (int i = 0; i < polynomialFunctions.length; i++){
            polynomialFunctions[i] = polynomialFunctions[i].multiply(constantPoly);
        }

        electronTemperature = new PolynomialSplineFunction(electronTemperature.getKnots(), polynomialFunctions);
    }


    /**
     * Simple toString() method useful for debugging
     */
    public String toString(){
        String string = " r (norm)  | rho (g/cc) |  Ti (keV)  |  Te (keV)\n";

        double[] r = getRadiusNodes();
        for (int i = 0; i < r.length; i++){
            string += String.format("%.4e | %.4e | %.4e | %.4e\n", r[i], massDensity.value(r[i]),
                    ionTemperature.value(r[i]), electronTemperature.value(r[i]));
        }

        return string;
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

        double getMagnitude() {
            return magnitude;
        }

        int getEl() {
            return el;
        }

        int getM() {
            return m;
        }

        void setMagnitude(double magnitude) {
            this.magnitude = magnitude;
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
     * Glorified struct to hold PlasmaLayer Species information
     */
    class PlasmaSpecies {

        private ParticleType type;
        private double numberFraction;

        PlasmaSpecies(int Z, double mass, double numberFraction) {
            this(new ParticleType(Z, mass), numberFraction);
        }

        public PlasmaSpecies(ParticleType type, double numberFraction) {
            this.type = type;
            this.numberFraction = numberFraction;
        }

        public ParticleType getType() {
            return type;
        }

        public double getZ(){
            return type.getZ();
        }

        public double getMass(){
            return type.getMass();
        }

        public double getNumberFraction() {
            return numberFraction;
        }
    }


    /**
     * Private convenience methods
     */

    private double getAverageMass(){
        double averageMass = 0.0;
        double totalFraction = 0.0;

        for (PlasmaSpecies species : plasmaSpecies){
            averageMass += species.getMass() * species.numberFraction;
            totalFraction += species.numberFraction;
        }

        return averageMass / totalFraction;
    }

    private double numberDensityFromRho(double rho){
        return rho / (getAverageMass() * 1000 * Constants.PROTON_MASS_KG);
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

    /**
     * Method that integrates over the volume of this PlasmaLayer weighting the results with weighting functions
     * that evaluate plasma properties at some position in the plasma
     */
    private double volumeIntegral(Method ... weightFunctions){
        try {
            final double MIN_PHI = 0.0;
            final double MAX_PHI = 2 * Math.PI;

            final double MIN_THETA = 0.0;
            final double MAX_THETA = 1 * Math.PI;

            double integral = 0.0;
            for (int k = 0; k < NUM_PHI_NODES; k++) {

                double phi_1 = MIN_PHI + (k + 1) * (MAX_PHI - MIN_PHI) / (NUM_PHI_NODES);
                double phi_0 = MIN_PHI + (k + 0) * (MAX_PHI - MIN_PHI) / (NUM_PHI_NODES);
                double dPhi = (phi_1 - phi_0);

                for (int j = 0; j < NUM_THETA_NODES; j++) {

                    double theta_1 = MIN_THETA + (j + 1) * (MAX_THETA - MIN_THETA) / (NUM_THETA_NODES);
                    double theta_0 = MIN_THETA + (j + 0) * (MAX_THETA - MIN_THETA) / (NUM_THETA_NODES);
                    double dTheta = (theta_1 - theta_0);

                    // Evaluate all of the outer radius bounds corresponding to this solid angle slice
                    double rMax_11 = getOuterRadiusBound(theta_1, phi_1);
                    double rMax_10 = getOuterRadiusBound(theta_1, phi_0);
                    double rMax_01 = getOuterRadiusBound(theta_0, phi_1);
                    double rMax_00 = getOuterRadiusBound(theta_0, phi_0);

                    // Evaluate all of the inner radius bounds corresponding to this solid angle slice
                    double rMin_11 = getInnerRadiusBound(theta_1, phi_1);
                    double rMin_10 = getInnerRadiusBound(theta_1, phi_0);
                    double rMin_01 = getInnerRadiusBound(theta_0, phi_1);
                    double rMin_00 = getInnerRadiusBound(theta_0, phi_0);


                    for (int i = 0; i < NUM_RADIUS_NODES; i++) {

                        // Determine the radius nodes corresponding to theta_1, phi_1
                        double r_111 = rMin_11 + (i + 1) * (rMax_11 - rMin_11) / NUM_RADIUS_NODES;
                        double r_011 = rMin_11 + (i + 0) * (rMax_11 - rMin_11) / NUM_RADIUS_NODES;
                        double dR_11 = (r_111 - r_011);

                        // Determine the radius nodes corresponding to theta_1, phi_0
                        double r_110 = rMin_10 + (i + 1) * (rMax_10 - rMin_10) / NUM_RADIUS_NODES;
                        double r_010 = rMin_10 + (i + 0) * (rMax_10 - rMin_10) / NUM_RADIUS_NODES;
                        double dR_10 = (r_110 - r_010);

                        // Determine the radius nodes corresponding to theta_0, phi_1
                        double r_101 = rMin_01 + (i + 1) * (rMax_01 - rMin_01) / NUM_RADIUS_NODES;
                        double r_001 = rMin_01 + (i + 0) * (rMax_01 - rMin_01) / NUM_RADIUS_NODES;
                        double dR_01 = (r_101 - r_001);

                        // Determine the radius nodes corresponding to theta_0, phi_0
                        double r_100 = rMin_00 + (i + 1) * (rMax_00 - rMin_00) / NUM_RADIUS_NODES;
                        double r_000 = rMin_00 + (i + 0) * (rMax_00 - rMin_00) / NUM_RADIUS_NODES;
                        double dR_00 = (r_100 - r_000);


                        // Create a vector for all 8 nodes
                        Vector3D vec_111 = Utils.getVectorFromSpherical(r_111, theta_1, phi_1);
                        Vector3D vec_011 = Utils.getVectorFromSpherical(r_011, theta_1, phi_1);

                        Vector3D vec_110 = Utils.getVectorFromSpherical(r_110, theta_1, phi_0);
                        Vector3D vec_010 = Utils.getVectorFromSpherical(r_010, theta_1, phi_0);

                        Vector3D vec_101 = Utils.getVectorFromSpherical(r_101, theta_0, phi_1);
                        Vector3D vec_001 = Utils.getVectorFromSpherical(r_001, theta_0, phi_1);

                        Vector3D vec_100 = Utils.getVectorFromSpherical(r_100, theta_0, phi_0);
                        Vector3D vec_000 = Utils.getVectorFromSpherical(r_000, theta_0, phi_0);


                        // Throw all the node quantities in arrays so we can loop through them easily
                        double[] rNodes = new double[] {r_111, r_011, r_110, r_010, r_101, r_001, r_100, r_000};
                        double[] dRs = new double[] {dR_11, dR_11, dR_10, dR_10, dR_01, dR_01, dR_00, dR_00};
                        double[] thetaNodes = new double[] {theta_1, theta_1, theta_1, theta_1, theta_0, theta_0, theta_0, theta_0};
                        Vector3D[] vectors = new Vector3D[] {vec_111, vec_011, vec_110, vec_010, vec_101, vec_001, vec_100, vec_000};


                        for (int nodeIndex = 0; nodeIndex < rNodes.length; nodeIndex++){

                            double value = rNodes[nodeIndex] * rNodes[nodeIndex];       // r^2 from Jacobian
                            value *= FastMath.sin(thetaNodes[nodeIndex]);               // sin(theta) from Jacobian
                            value *= dRs[nodeIndex] * dTheta * dPhi / 8.0;              // Integration widths

                            // Multiply all of the weight functions
                            for (Method weight : weightFunctions){
                                value *= (double) weight.invoke(this, vectors[nodeIndex]);
                            }

                            integral += value;
                        }
                    }
                }
            }

            return integral;

        }catch (Exception e){
            e.printStackTrace();
            System.exit(-1);
        }

        return Double.NaN;
    }

}