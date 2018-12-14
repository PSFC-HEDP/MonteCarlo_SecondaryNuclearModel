package MonteCarloParticleTracer;

import cStopPow.DoubleVector;
import cStopPow.StopPow_LP;
import cStopPow.StopPow_Zimmerman;
import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Method;
import java.util.ArrayList;

/**
 * This class represents a single plasma of plasma defined by radial profiles of mass density, ion temperature, and
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
public class Plasma {

    // Constants used for volumetric integration
    // Benchmark note (50/50/50 gets sphere volume correct with 0.013%)
    private static final int NUM_PHI_NODES    = 50;
    private static final int NUM_THETA_NODES  = 50;
    private static final int NUM_RADIUS_NODES = 50;

    // Constants for handling floating point exceptions
    private static final double OUT_BOUNDS_EPSILON = 1e-6;


    // ArrayList of Legendre Modes that describe the inner boundary of this Plasma
    private ArrayList<LegendreMode> innerBoundaryLegendreModes = new ArrayList<>();

    // ArrayList of Legendre Modes that describe the outer boundary of this Plasma
    private ArrayList<LegendreMode> outerBoundaryLegendreModes = new ArrayList<>();

    // ArrayList of Particle Species in this Plasma
    private ArrayList<PlasmaSpecies> plasmaSpecies = new ArrayList<>();

    // Ion Temperature radial profile define between [rMin, rMax] (keV)
    // Internally, binEdges are mapped to the normalized units of [0, 1]
    private PolynomialSplineFunction ionTemperature;

    // Electron Temperature radial profile defined between [rMin, rMax] (keV)
    // Internally, binEdges are mapped to the normalized units of [0, 1]
    private PolynomialSplineFunction electronTemperature;

    // Mass density radial profile defined between [rMin, rMax] (keV)
    // Internally, binEdges are mapped to the normalized units of [0, 1]
    private PolynomialSplineFunction massDensity;


    // Some flags we'll use for special cases that can be used to speed up computational time
    private boolean boundsVaryInTheta = false;
    private boolean boundsVaryInPhi   = false;



    /**
     * Default constructor for building a plasma object
     * At least 1 species needs to be added in addition to this construction
     * At least 1 outer Legendre mode (usually a P0) needs to be defined in addition to this construction
     *
     * @param radiusNodes Nodes (units arbitrary) on which the plasma properties are known (will be normalized internally)
     * @param ionTemperature Ion temperature (in keV) evaluated at the radius nodes
     * @param electronTemperature Electron temperature (in keV) evaluated at the radius nodes
     * @param massDensity Mass density (in g/cc) evaluated at the radius nodes
     */
    public Plasma(double[] radiusNodes, double[] ionTemperature, double[] electronTemperature, double[] massDensity) {
        // Normalize the radii to [0, 1]
        DoubleArray normalizedRadii = new DoubleArray(radiusNodes);
        normalizedRadii.subtract(normalizedRadii.getMin());
        normalizedRadii.multiply(1.0 / normalizedRadii.getMax());

        // Floating point issues, we need to edges to be EXACTLY [0, 1]
        normalizedRadii.set(0, 0);
        normalizedRadii.set(normalizedRadii.length()-1, 1);

        this.ionTemperature = new SplineInterpolator().interpolate(normalizedRadii.getValues(), ionTemperature);
        this.electronTemperature = new SplineInterpolator().interpolate(normalizedRadii.getValues(), electronTemperature);
        this.massDensity = new SplineInterpolator().interpolate(normalizedRadii.getValues(), massDensity);
    }


    /**
     * Returns a plasma with a constant Ti = Te = 1 keV and constant rho = 1 g/cc
     * Use setBurnAveragedIonTemperature to set ion temperature
     * Use setElectronTemperatureFraction to set electron temperature
     * Use setAverageMassDensity if mass density is known
     * Use setTotalMass if total mass is known
     *
     * @param numberOfNodes nodes for the profiles
     * @return a uniform plasma
     */
    public static Plasma uniformPlasma(int numberOfNodes){
        double[] r = DoubleArray.linspace(0, 1, numberOfNodes).getValues();
        double[] ones = DoubleArray.linspace(1, 1, numberOfNodes).getValues();
        return new Plasma(r, ones, ones, ones);
    }

    public Plasma copy(){

        // Recreate the profiles
        double[] r   = this.getNormalizedRadiusNodes();
        double[] Ti  = new double[r.length];
        double[] Te  = new double[r.length];
        double[] rho = new double[r.length];
        for (int i = 0; i < r.length; i++){
            Ti[i]  = ionTemperature.value(r[i]);
            Te[i]  = electronTemperature.value(r[i]);
            rho[i] = massDensity.value(r[i]);
        }

        // Create the plasma object
        Plasma copy = new Plasma(r, Ti, Te, rho);

        // Inner Legendre Modes
        for (LegendreMode mode : innerBoundaryLegendreModes){
            copy.addInnerLegendreMode(mode.getEl(), mode.getM(), mode.getMagnitude());
        }

        // Outer Legendre Modes
        for (LegendreMode mode : outerBoundaryLegendreModes){
            copy.addOuterLegendreMode(mode.getEl(), mode.getM(), mode.getMagnitude());
        }

        // Plasma Species
        for (PlasmaSpecies species : plasmaSpecies){
            copy.addSpecies(new ParticleType(species.getZ(), species.getMass()),
                    species.getNumberFraction());
        }

        return copy;
    }



    // ***************************************************
    // Add methods used to aid in constructing the Plasma
    // ***************************************************


    public void addInnerLegendreMode(int l, int m, double magnitude){

        // Update our flags if needed
        if (m != 0) boundsVaryInPhi = true;
        if (l != 0) boundsVaryInTheta = true;

        // Add this mode
        innerBoundaryLegendreModes.add(new LegendreMode(l, m, magnitude));

    }

    public void addOuterLegendreMode(int l, int m, double magnitude){

        // Update our flags if needed
        if (m != 0) boundsVaryInPhi = true;
        if (l != 0) boundsVaryInTheta = true;

        // Add this mode
        outerBoundaryLegendreModes.add(new LegendreMode(l, m, magnitude));

    }

    public void addSpecies(ParticleType type, double numberProportion){
        plasmaSpecies.add(new PlasmaSpecies(type, numberProportion));

        // Re-normalize the number fractions
        double totalProportion = 0.0;
        for (PlasmaSpecies species : plasmaSpecies){
            totalProportion += species.numberProportion;
        }

        for (PlasmaSpecies species : plasmaSpecies){
            species.numberFraction = species.numberProportion / totalProportion;
        }
    }

    // **************
    // Getter methods
    // **************

    /**
     * Method that returns an array of values corresponding the the stopping power (dE/dx) of the testParticle
     * at every radius normalized location in this Plasma in units of MeV / cm
     * @param testParticle ParticleType whose dEdx is to be calculated (depends on Z and A)
     * @param energy Energy (in MeV) of the particle whose dEdx is to be calculated
     * @return An array of dEdx(r) with the same length as the radial nodes used to define the profiles
     */
    double[] getStoppingPower(ParticleType testParticle, double energy){

        double[] r = getNormalizedRadiusNodes();                  // Normalized radius [0, 1]
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



            if (testParticle.getZ() == 0){
                dEdx[i] = 0.0;
            }else {

                /*
                StopPow_Zimmerman stopPow_zimmerman = new StopPow_Zimmerman(testParticle.getMass(), testParticle.getZ(),
                        speciesAs, speciesZs, speciesTs, speciesNs, speciesZbars, Te);
                dEdx[i] = 1e4 * stopPow_zimmerman.dEdx_MeV_um(energy);
                */

                StopPow_LP stopPow_lp = new StopPow_LP(testParticle.getMass(), testParticle.getZ(),
                        speciesAs, speciesZs, speciesTs, speciesNs, Te);
                dEdx[i] = 1e4 * stopPow_lp.dEdx_MeV_um(energy);

            }
        }

        return dEdx;
    }

    double[] getNormalizedRadiusNodes(){
        return massDensity.getKnots();
    }

    ArrayList<LegendreMode> getInnerBoundaryLegendreModes() {
        return innerBoundaryLegendreModes;
    }

    ArrayList<LegendreMode> getOuterBoundaryLegendreModes() {
        return outerBoundaryLegendreModes;
    }

    boolean boundsVaryInTheta(){
        return boundsVaryInTheta;
    }

    boolean boundsVaryInPhi(){
        return boundsVaryInPhi;
    }



    // ******************************************************
    // Getter methods that vary as a function of (theta, phi)
    // ******************************************************

    double getInnerRadiusBound(double theta, double phi){
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

    public void DEBUG_dumpOuterBounds(){
        try {
            FileWriter w = new FileWriter("bounds.csv");

            double[] theta = DoubleArray.linspace(0, Math.PI, 501).getValues();
            double[] phi   = DoubleArray.linspace(0, 2*Math.PI, 501).getValues();

            for (int i = 0; i < theta.length; i++){
                w.write("," + theta[i]);
            }
            w.write("\n");

            for (int j = 0; j < phi.length; j++){
                w.write("" + phi[j]);

                for (int i = 0; i < theta.length; i++){
                    w.write("," + getOuterRadiusBound(theta[i], phi[i]));
                }
                w.write("\n");
            }

            w.close();

        }
        catch (IOException e){
            e.printStackTrace();
        }

    }

    public void DEBUG_printOuterModes(){
        for (LegendreMode mode : outerBoundaryLegendreModes){
            System.out.println(mode.getEl() + " " + mode.getM() + " " + mode.getMagnitude() + " " + mode.getNorm() * mode.getMagnitude());
        }
    }

    double[] getRadiusNodes(double theta, double phi){
        double[] rN = getNormalizedRadiusNodes();
        double rMin = getInnerRadiusBound(theta, phi);
        double rMax = getOuterRadiusBound(theta, phi);

        double[] r = new double[rN.length];
        for (int i = 0; i < rN.length; i++){
            r[i] = rMin + (rMax - rMin)*rN[i];
        }

        return r;
    }

    public double getArealDensity(double theta, double phi){

        TrapezoidIntegrator integrator = new TrapezoidIntegrator();
        double value = integrator.integrate(Integer.MAX_VALUE, massDensity, 0, 1);

        double Rmin = getInnerRadiusBound(theta, phi);
        double Rmax = getOuterRadiusBound(theta, phi);

        // We have to multiply by this normalization factor since the profiles are defined in normalized radius
        return (Rmax - Rmin) * value;
    }


    // ***********************************************************************************
    // Getter methods that evaluate the plasma conditions at some specific position vector
    // ***********************************************************************************

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

            e.printStackTrace();
            System.err.printf("Warning! Ion temperature evaluation %.4f%% out of bounds\n", 100*(R - maxRadius)/maxRadius);
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

            System.err.printf("Warning! Electron temperature evaluation %.4f%% out of bounds\n", 100*(R - maxRadius)/maxRadius);
        }

        return Double.NaN;
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

            System.err.printf("Warning! Mass density evaluation %.4f%% out of bounds\n", 100*(R - maxRadius)/maxRadius);
        }

        return Double.NaN;
    }

    public double getSpeciesMassDensity(Vector3D r, ParticleType type){
        return getSpeciesNumberDensity(r, type) * type.getMass() * Constants.GRAMS_PER_AMU;
    }

    public double getIonNumberDensity(Vector3D r){
        return numberDensityFromRho(getMassDensity(r));
    }

    public double getElectronNumberDensity(Vector3D r){

        // Init the electron density
        double ne = 0;

        // Get the total ion density
        double n   = getIonNumberDensity(r);

        // Loop through all the plasma species
        for (PlasmaSpecies species : plasmaSpecies){

            // Get the number fraction
            double fi = species.getNumberFraction();

            // Get the Zbar (using Drake ionization assumption)
            double Zbar = 20*Math.sqrt(getIonTemperature(r));
            Zbar = Math.min(Zbar, species.getZ());

            // Add these electrons to the total
            ne += (fi * n * Zbar);
        }

        // Return
        return ne;

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

        return getIonNumberDensity(r) * speciesFraction / totalFraction;
    }

    double getThermalSigma(Vector3D position, NuclearReaction reaction){
        // TODO: Assuming 2 body
        ParticleType[] reactants = reaction.getReactants();
        ParticleType A = reactants[0];
        ParticleType B = reactants[1];

        ParticleType[] products = reaction.getProducts();
        ParticleType C = products[0];
        ParticleType D = products[1];

        // Calculate factor = 2 * mC * mD / (mA + mB) / (mC + mD)
        double factor = 2.0 * C.getMass() * D.getMass();
        factor /= (A.getMass() + B.getMass());
        factor /= (C.getMass() + D.getMass());

        // Get the Q value
        double Q = reaction.getQValue();

        // Get the ion temperature
        double Tion_MeV = 1e-3 * getIonTemperature(position);

        // Return the sigma in MeV
        return Math.sqrt(factor * Q * Tion_MeV);
    }

    boolean containsSpecies(ParticleType type){

        for (PlasmaSpecies species : plasmaSpecies){
            if (species.getType().equals(type))  return true;
        }

        return false;
    }



    // ********************************************************************
    // Getter methods for evaluating the reactivity at some position vector
    // ********************************************************************

    public double getReactivity(Vector3D r, Reactivity reactivity){
        return reactivity.evaluate(getIonTemperature(r));
    }

    public double getDDnReactivity(Vector3D r){
        return getReactivity(r, Reactivity.DDn_Reactivity);
    }

    public double getDDpReactivity(Vector3D r){
        return getReactivity(r, Reactivity.DDp_Reactivity);
    }

    public double getD3HepReactivity(Vector3D r){
        return getReactivity(r, Reactivity.D3Hep_Reactivity);
    }

    public double get3He3HepReactivity(Vector3D r){ return getReactivity(r, Reactivity.He3He3_Reactivity); }



    // **************************************************************
    // Methods used for handling particles near the plasma's boundary
    // **************************************************************

    boolean getIsInside(Vector3D position){
        double[] coordinates = Utils.getSphericalFromVector(position);
        double theta = coordinates[1];
        double phi   = coordinates[2];

        if (position.getNorm() < getOuterRadiusBound(theta, phi)){

            // This plasma may not have an inner bound
            if (innerBoundaryLegendreModes.size() == 0)
                return true;

            // Otherwise, verify that we're within the inner bound as well
            if (position.getNorm() > getInnerRadiusBound(theta, phi))
                return true;

        }

        return false;
    }

    Particle moveInside(Particle particle){
        if (this.getIsInside(particle.getPosition())){
            return particle;
        }

        double[] coordinates = Utils.getSphericalFromVector(particle.getPosition());
        double theta = coordinates[1];
        double phi   = coordinates[2];

        double radius;
        if (this.getIsCloserToOuterBoundary(particle.getPosition())){
            radius = getOuterRadiusBound(theta, phi);
        }else{
            radius = getInnerRadiusBound(theta, phi);
        }

        particle.setPosition(Utils.getVectorFromSpherical(radius, theta, phi));
        return particle;
    }

    boolean getIsCloserToOuterBoundary(Vector3D r){
        double[] coordinates = Utils.getSphericalFromVector(r);
        double theta = coordinates[1];
        double phi   = coordinates[2];

        double distanceToOuter = getOuterRadiusBound(theta, phi) - r.getNorm();
        double distanceToInner = r.getNorm() - getInnerRadiusBound(theta, phi);

        return distanceToOuter < distanceToInner;
    }

    boolean getIsCloserToInnerBoundary(Vector3D r){
        return !getIsCloserToOuterBoundary(r);
    }



    // ****************************************
    // Getters for volume integrated quantities
    // ****************************************

    double getTotalVolume(){
        return volumeIntegral(new String[] {});
    }

    public double getTotalMass(){
        return volumeIntegral("getMassDensity");
    }

    double getVolumeAveragedIonTemperature(){
        return volumeIntegral("getIonTemperature") / getTotalVolume();
    }

    double getVolumeAveragedElectronTemperature(){
        return volumeIntegral("getElectronTemperature") / getTotalVolume();
    }

    double getVolumeAveragedMassDensity(){
        return volumeIntegral("getMassDensity") / getTotalVolume();
    }



    // ************************
    // Burn rate getter methods
    // ************************

    private double getBurnRate(String reactivityMethodName){
        return volumeIntegral("getMassDensity", "getMassDensity", reactivityMethodName);
    }

    public double getDDnBurnRate(){
        return getBurnRate("getDDnReactivity");
    }

    public double getDDpBurnRate(){
        return getBurnRate("getDDpReactivity");
    }

    public double getD3HepBurnRate(){
        return getBurnRate("getD3HepReactivity");
    }

    public double get3He3HeBurnRate(){
        return getBurnRate("get3He3HepReactivity");
    }



    // ****************************
    // Burn averaged getter methods
    // ****************************

    private double getBurnAveragedIonTemperature(String reactivityMethodName){
        double numerator   = volumeIntegral("getMassDensity", "getMassDensity", reactivityMethodName, "getIonTemperature");
        double denominator = volumeIntegral("getMassDensity", "getMassDensity", reactivityMethodName);
        return  numerator / denominator;
    }

    public double getDDnBurnAveragedIonTemperature(){
        return getBurnAveragedIonTemperature("getDDnReactivity");
    }

    public double getDDpBurnAveragedIonTemperature(){
        return getBurnAveragedIonTemperature("getDDpReactivity");
    }

    public double getD3HepBurnAveragedIonTemperature(){
        return getBurnAveragedIonTemperature("getD3HepReactivity");
    }

    public double get3He3HepBurnAveragedIonTemperature(){ return getBurnAveragedIonTemperature("get3He3HepReactivity"); }


    private double getBurnAveragedElectronTemperature(String reactivityMethodName){
        double numerator   = volumeIntegral("getMassDensity", "getMassDensity", reactivityMethodName, "getElectronTemperature");
        double denominator = volumeIntegral("getMassDensity", "getMassDensity", reactivityMethodName);
        return  numerator / denominator;
    }

    public double getDDnBurnAveragedElectronTemperature(){
        return getBurnAveragedElectronTemperature("getDDnReactivity");
    }

    public double getDDpBurnAveragedElectronTemperature(){
        return getBurnAveragedElectronTemperature("getDDpReactivity");
    }

    public double getD3HepBurnAveragedElectronTemperature(){
        return getBurnAveragedElectronTemperature("getD3HepReactivity");
    }

    public double get3He3HepBurnAveragedElectronTemperature(){
        return getBurnAveragedElectronTemperature("get3He3HepReactivity");
    }


    // ********************************
    // Burn distribution getter methods
    // ********************************

    Distribution getRadialBremDistribution(){

        double[] r = getNormalizedRadiusNodes();                // Normalized r
        double[] J = new double[r.length];                      // Brem probability

        for (int i = 0; i < r.length; i++) {

            // Get the plasma conditions at this position
            double n = numberDensityFromRho(massDensity.value(r[i]));
            double Ti = ionTemperature.value(r[i]);
            double Te = electronTemperature.value(r[i]);


            // Loop through the plasma species
            double ni_Z_2 = 0, ne = 0;
            for (PlasmaSpecies species : plasmaSpecies) {

                // Get the number fraction
                double fi = species.getNumberFraction();

                // Get the Zbar (using Drake ionization assumption)
                double Zbar = 20 * Math.sqrt(Ti);
                Zbar = Math.min(Zbar, species.getZ());

                // Add this species to the n_i Z^2 term
                ni_Z_2 += (fi * n * Zbar * Zbar);

                // Add this species to the n_i Z (ne) term
                ne += (fi * n * Zbar);

            }

            // Brem is (ne ni Z^2) * Te^(1/2). r^2 factor because spherical
            J[i] = r[i] * r[i] * ne * ni_Z_2 * Math.sqrt(Te);

        }

        return new Distribution(r, J);

    }

    Distribution getRadialBurnDistribution(Reactivity reactivity){

        double[] r = getNormalizedRadiusNodes();                  // Normalized r
        double[] Y = new double[r.length];              // Yield probability

        for (int i = 0; i < r.length; i++){
            double rho = massDensity.value(r[i]);
            double T   = ionTemperature.value(r[i]);
            double sigmaV = reactivity.evaluate(T);

            Y[i] = r[i] * r[i] * rho * rho * sigmaV;        // Y(r) \proto n^2 <\sigma v>. r^2 factor because spherical
        }

        return new Distribution(r, Y);
    }

    Distribution getRadialDDnBurnDistribution(){
        return getRadialBurnDistribution(Reactivity.DDn_Reactivity);
    }

    Distribution getRadialDDpBurnDistribution(){
        return getRadialBurnDistribution(Reactivity.DDp_Reactivity);
    }

    Distribution getRadialD3HepBurnDistribution(){
        return getRadialBurnDistribution(Reactivity.D3Hep_Reactivity);
    }

    Distribution getRadial3He3HepBurnDistribution(){
        return getRadialBurnDistribution(Reactivity.He3He3_Reactivity);
    }

    Distribution getPolarDistribution(){

        final int NUM_THETA_NODES = 201;
        double[] theta  = DoubleArray.linspace(1e-6, Math.PI-1e-6, NUM_THETA_NODES).getValues();
        double[] probability = new double[theta.length];

        for (int i = 0; i < theta.length; i++){
            double Rmin = getInnerRadiusBound(theta[i], 0.0);
            double Rmax = getOuterRadiusBound(theta[i], 0.0);
            probability[i] = Math.sin(theta[i]) * (Rmax - Rmin);
        }

        return new Distribution(theta, probability);

    }




    // **************
    // Setter methods
    // **************

    void setInnerBoundaryLegendreModes(ArrayList<LegendreMode> innerBoundaryLegendreModes) {
        this.innerBoundaryLegendreModes = new ArrayList<>();
        for (LegendreMode mode : innerBoundaryLegendreModes){
            this.addInnerLegendreMode(mode.getEl(), mode.getM(), mode.getMagnitude());
        }
    }

    void setOuterBoundaryLegendreModes(ArrayList<LegendreMode> outerBoundaryLegendreModes) {
        this.outerBoundaryLegendreModes = new ArrayList<>();
        for (LegendreMode mode : outerBoundaryLegendreModes){
            this.addInnerLegendreMode(mode.getEl(), mode.getM(), mode.getMagnitude());
        }
    }

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

    public void setElectronTemperatureFraction(double electronTemperatureFraction){
        double currentFraction = getVolumeAveragedElectronTemperature() / getVolumeAveragedIonTemperature();
        double multiplicativeFactor = electronTemperatureFraction / currentFraction;
        multiplyElectronTemperatureByScalar(multiplicativeFactor);
    }

    public void setTotalMass(double totalMass) {
        double multiplicativeFactor = totalMass / getTotalMass();
        multiplyMassDensityByScalar(multiplicativeFactor);
    }

    public void setAverageMassDensity(double massDensity){
        double multiplicativeFactor = massDensity / getVolumeAveragedMassDensity();
        multiplyMassDensityByScalar(multiplicativeFactor);
    }



    // ****************************
    // Set burn averaged quantities
    // ****************************

    private void setBurnAveragedIonTemperature(double burnT, String reactivityMethodName){
        double multiplicativeFactor = burnT / getBurnAveragedIonTemperature(reactivityMethodName);
        multiplyIonTemperatureByScalar(multiplicativeFactor);

    }

    public void setDDnBurnAveragedIonTemperature(double burnT){
        setBurnAveragedIonTemperature(burnT, "getDDnReactivity");
    }

    public void setDDpBurnAveragedIonTemperature(double burnT){
        setBurnAveragedIonTemperature(burnT, "getDDpReactivity");
    }

    public void setD3HepBurnAveragedIonTemperature(double burnT){
        setBurnAveragedIonTemperature(burnT, "getD3HepReactivity");
    }

    public void set3He3HeBurnAveragedIonTemperature(double burnT){
        setBurnAveragedIonTemperature(burnT, "get3He3HepReactivity");
    }



    // *******************
    // toString() Override
    // *******************

    public String toString(){
        String string = "r (norm), rho (g/cc), Ti (keV), Te (keV)\n";

        double[] r = getNormalizedRadiusNodes();
        for (int i = 0; i < r.length; i++){
            string += String.format("%.8e, %.8e, %.8e, %.8e\n", r[i], massDensity.value(r[i]),
                    ionTemperature.value(r[i]), electronTemperature.value(r[i]));
        }

        return string;
    }



    /**
     * Internal Legendre Mode class
     */



    /**
     * Internal plasma species class
     * Glorified struct that just pairs a ParticleType and a number fraction
     */
    private class PlasmaSpecies {

        private ParticleType type;
        private double numberProportion;    // Not guaranteed to be normalized
        private double numberFraction;      // Guaranteed to be normalized

        PlasmaSpecies(ParticleType type, double numberProportion) {
            this.type = type;
            this.numberProportion = numberProportion;
            this.numberFraction = numberProportion;
        }

        ParticleType getType() {
            return type;
        }

        int getZ(){
            return type.getZ();
        }

        double getMass(){
            return type.getMass();
        }

        double getNumberFraction() {
            return numberFraction;
        }
    }



    // ****************************************************************
    // Private methods for multiplying the profiles by scalar constants
    // ****************************************************************

    public void multiplyMassDensityByScalar(double scalar){
        PolynomialFunction constantPoly = new PolynomialFunction(new double[] {scalar});

        PolynomialFunction[] polynomialFunctions = massDensity.getPolynomials();
        for (int i = 0; i < polynomialFunctions.length; i++){
            polynomialFunctions[i] = polynomialFunctions[i].multiply(constantPoly);
        }

        massDensity = new PolynomialSplineFunction(massDensity.getKnots(), polynomialFunctions);
    }

    public void multiplyIonTemperatureByScalar(double scalar){
        PolynomialFunction constantPoly = new PolynomialFunction(new double[] {scalar});

        PolynomialFunction[] polynomialFunctions = ionTemperature.getPolynomials();
        for (int i = 0; i < polynomialFunctions.length; i++){
            polynomialFunctions[i] = polynomialFunctions[i].multiply(constantPoly);
        }

        ionTemperature = new PolynomialSplineFunction(ionTemperature.getKnots(), polynomialFunctions);
    }

    public void multiplyElectronTemperatureByScalar(double scalar){
        PolynomialFunction constantPoly = new PolynomialFunction(new double[] {scalar});

        PolynomialFunction[] polynomialFunctions = electronTemperature.getPolynomials();
        for (int i = 0; i < polynomialFunctions.length; i++){
            polynomialFunctions[i] = polynomialFunctions[i].multiply(constantPoly);
        }

        electronTemperature = new PolynomialSplineFunction(electronTemperature.getKnots(), polynomialFunctions);
    }



    // *********************************************
    // Various other private methods for convenience
    // *********************************************

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



    // ***********************
    // Volume integral methods
    // ***********************

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

                        // Determine the radius binEdges corresponding to theta_1, phi_1
                        double r_111 = rMin_11 + (i + 1) * (rMax_11 - rMin_11) / NUM_RADIUS_NODES;
                        double r_011 = rMin_11 + (i + 0) * (rMax_11 - rMin_11) / NUM_RADIUS_NODES;
                        double dR_11 = (r_111 - r_011);

                        // Determine the radius binEdges corresponding to theta_1, phi_0
                        double r_110 = rMin_10 + (i + 1) * (rMax_10 - rMin_10) / NUM_RADIUS_NODES;
                        double r_010 = rMin_10 + (i + 0) * (rMax_10 - rMin_10) / NUM_RADIUS_NODES;
                        double dR_10 = (r_110 - r_010);

                        // Determine the radius binEdges corresponding to theta_0, phi_1
                        double r_101 = rMin_01 + (i + 1) * (rMax_01 - rMin_01) / NUM_RADIUS_NODES;
                        double r_001 = rMin_01 + (i + 0) * (rMax_01 - rMin_01) / NUM_RADIUS_NODES;
                        double dR_01 = (r_101 - r_001);

                        // Determine the radius binEdges corresponding to theta_0, phi_0
                        double r_100 = rMin_00 + (i + 1) * (rMax_00 - rMin_00) / NUM_RADIUS_NODES;
                        double r_000 = rMin_00 + (i + 0) * (rMax_00 - rMin_00) / NUM_RADIUS_NODES;
                        double dR_00 = (r_100 - r_000);


                        // Create a vector for all 8 binEdges
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