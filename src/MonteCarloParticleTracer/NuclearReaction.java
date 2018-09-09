package MonteCarloParticleTracer;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.FastMath;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * Class that handles a 2 body nuclear reaction A + B -> C + D
 */
public class NuclearReaction {

    private ParticleType[] reactants;       // Reacting particles
    private ParticleType[] products;        // Resulting products

    private PolynomialSplineFunction crossSection;      // Function representing the XS in cm^2 as a function of CoM energy in MeV

    public static final NuclearReaction DD_t   = new NuclearReaction(ParticleType.deuteron, ParticleType.deuteron, ParticleType.triton, ParticleType.proton);
    public static final NuclearReaction DD_3He = new NuclearReaction(ParticleType.deuteron, ParticleType.deuteron, ParticleType.helium3, ParticleType.neutron);

    public static final NuclearReaction DT_n = new NuclearReaction(ParticleType.deuteron, ParticleType.triton, ParticleType.neutron, ParticleType.alpha, DataFiles.DT_XS_File);
    public static final NuclearReaction D3He_p = new NuclearReaction(ParticleType.deuteron, ParticleType.helium3, ParticleType.proton, ParticleType.alpha, DataFiles.D3He_XS_File);

    public NuclearReaction(ParticleType A, ParticleType B, ParticleType C, ParticleType D, File crossSectionFile) {

        // Set the reactants
        this.reactants = new ParticleType[2];
        this.reactants[0] = A;
        this.reactants[1] = B;

        // Set the products
        this.products = new ParticleType[2];
        this.products[0] = C;
        this.products[1] = D;

        this.crossSection = generateCrossSectionFunction(crossSectionFile);
    }

    public NuclearReaction(ParticleType A, ParticleType B, ParticleType C, ParticleType D) {

        // Set the reactants
        this.reactants = new ParticleType[2];
        this.reactants[0] = A;
        this.reactants[1] = B;

        // Set the products
        this.products = new ParticleType[2];
        this.products[0] = C;
        this.products[1] = D;

    }

    public NuclearReaction copy(){

        // Handle the reactants
        ParticleType A = new ParticleType(this.reactants[0].getZ(), this.reactants[0].getMass());
        ParticleType B = new ParticleType(this.reactants[1].getZ(), this.reactants[1].getMass());

        // Handled the products
        ParticleType C = new ParticleType(this.products[0].getZ(), this.products[0].getMass());
        ParticleType D = new ParticleType(this.products[1].getZ(), this.products[1].getMass());

        // Create the copy
        NuclearReaction copy = new NuclearReaction(A, B, C, D);

        // If the cross section exists, add that
        if (crossSection != null){

            double[] energies = crossSection.getKnots();
            double[] values   = new double[energies.length];
            for (int i = 0; i < energies.length; i++){
                values[i] = crossSection.value(energies[i]);
            }

            copy.crossSection = new SplineInterpolator().interpolate(energies, values);
        }

        return copy;
    }

    // TODO: This assumes 2 body
    public Particle getProductParticle(Particle A, Particle B, Vector3D direction){

        if (direction == null){
            return getProductParticle(A, B);
        }

        // Rename products for readability
        ParticleType C = this.products[0];      // The returned product
        ParticleType D = this.products[1];      // The other product

        // Reduced mass of these two particles in amu
        double reducedMass = A.getMass()*B.getMass() / (A.getMass() + B.getMass());

        // Relative velocity of these two particles as a fraction of c
        Vector3D relativeVelocity = A.getVelocity().subtract(B.getVelocity());

        // Center of mass kinetic energy of this system in MeV
        double kineticEnergy = 0.5*Constants.MEV_PER_AMU*reducedMass*FastMath.pow(relativeVelocity.getNorm(),2);    // MeV

        // Center of mass velocity of this system in MeV
        Vector3D centerOfMassVelocity = getCenterOfMassVelocity(A, B);

        // Center of mass energy of the product of interest in MeV
        double productEnergy_CoM = (D.getMass() / (C.getMass() + D.getMass())) * (getQValue() + kineticEnergy);

        // Center of mass speed of the product of interest as a fraction of c
        double productSpeed_CoM = FastMath.sqrt(2*productEnergy_CoM/C.getMass()/Constants.MEV_PER_AMU);


        /**
         * How we handle it when we DO force the direction
         */

        // Angle between the center of mass velocity and the forced direction
        double angle = Vector3D.angle(centerOfMassVelocity, direction);

        // Use the law of cosines (based on vC = vCM + uC) to get the speed of the product as a fraction of c
        double productSpeed_Lab = FastMath.pow(centerOfMassVelocity.getNorm(), 2);       // vCM^2
        productSpeed_Lab *= FastMath.pow(FastMath.cos(angle),2) - 1;                     // vCM^2 * (cos^2(\theta) - 1)
        productSpeed_Lab += FastMath.pow(productSpeed_CoM,2);                            // uC^2 + vCM^2 * (cos^2(\theta) - 1)
        productSpeed_Lab = FastMath.sqrt(productSpeed_Lab);                                 // \sqrt(uC^2 + vCM^2 * (cos^2(\theta) - 1))
        productSpeed_Lab += centerOfMassVelocity.getNorm()*FastMath.cos(angle);             // vCM * cos(\theta) + \sqrt(uC^2 + vCM^2 * (cos^2(\theta) - 1))

        // Convert this to an energy in MeV
        double productEnergy_Lab = 0.5*Constants.MEV_PER_AMU*C.getMass()*FastMath.pow(productSpeed_Lab,2);        // MeV


        // Build the  particle of type C born at the same position as A with the forced direction and our derived energy
        Particle productParticle = new Particle(C, A.getPosition(), direction, productEnergy_Lab, A.getTime());

        // Normalization associated with forcing direction
        productParticle.multiplyWeight(productEnergy_Lab / productEnergy_CoM);

        return productParticle;
    }

    // TODO: This assumes 2 body
    private Particle getProductParticle(Particle A, Particle B){

        // Rename products for readability
        ParticleType C = this.products[0];      // The returned product
        ParticleType D = this.products[1];      // The other product

        // Reduced mass of these two particles in amu
        double reducedMass = A.getMass()*B.getMass() / (A.getMass() + B.getMass());

        // Relative velocity of these two particles as a fraction of c
        Vector3D relativeVelocity = A.getVelocity().subtract(B.getVelocity());

        // Center of mass kinetic energy of this system in MeV
        double kineticEnergy = 0.5*Constants.MEV_PER_AMU*reducedMass*FastMath.pow(relativeVelocity.getNorm(),2);    // MeV

        // Center of mass velocity of this system in MeV
        Vector3D centerOfMassVelocity = getCenterOfMassVelocity(A, B);

        // Center of mass energy of the product of interest in MeV
        double productEnergy_CoM = (D.getMass() / (C.getMass() + D.getMass())) * (getQValue() + kineticEnergy);

        // Center of mass speed of the product of interest as a fraction of c
        double productSpeed_CoM = FastMath.sqrt(2*productEnergy_CoM/C.getMass()/Constants.MEV_PER_AMU);


        /**
         * How we handle it if we didn't force the direction
         */

        Vector3D productVelocity_CoM = Utils.sampleRandomNormalizedVector().scalarMultiply(productSpeed_CoM);
        Vector3D productVelocity_Lab = centerOfMassVelocity.add(productVelocity_CoM);
        double productSpeed_Lab = productVelocity_Lab.getNorm();
        double productEnergy_Lab = 0.5*Constants.MEV_PER_AMU*C.getMass()*FastMath.pow(productSpeed_Lab,2);        // MeV

        Particle productParticle = new Particle(C, A.getPosition(), productVelocity_Lab.normalize(), productEnergy_Lab, A.getTime());


        return productParticle;
    }



    // Getter functions
    // ****************

    /**
     * Calculate and return the Q of this reaction
     * @return Q of this reaction in MeV
     */
    public Double getQValue(){

        // Initialize our Q
        double Q = 0.0;

        // Add the mass of all the reacting particles
        for (ParticleType reactant : reactants){
            Q += reactant.getMass();
        }

        // Subtract the mass of all of our product particles
        for (ParticleType product : products){
            Q -= product.getMass();
        }

        // Return in MeV
        return Q * Constants.MEV_PER_AMU;
    }

    /**
     * Returns the maximum possible energy product
     * @param A
     * @param B
     * @return
     */
    public Particle getMaxEnergyProductParticle(Particle A, Particle B){
        Vector3D centerOfMassVelocity = getCenterOfMassVelocity(A, B);
        Vector3D direction = centerOfMassVelocity.normalize();
        return getProductParticle(A, B, direction);
    }

    public Particle getMinEnergyProductParticle(Particle A, Particle B){
        Vector3D centerOfMassVelocity = getCenterOfMassVelocity(A, B);
        Vector3D direction = centerOfMassVelocity.scalarMultiply(-1.0).normalize();
        return getProductParticle(A, B, direction);
    }

    public double getZeroTemperatureMeanEnergy(){
        ParticleType C = products[0];
        ParticleType D = products[1];

        return getQValue() * D.getMass() / (C.getMass() + D.getMass());
    }

    public ParticleType[] getReactants() {
        return reactants;
    }

    public ParticleType[] getProducts() {
        return products;
    }

    public String toString(){
        String string = "";

        // Add all of the reactants
        for (ParticleType reactant : reactants) {
            string += reactant.toString() + " + ";
        }

        // Remove the extra " + "
        string = string.substring(0, string.length() - 3);

        // Add the " -> "
        string += " -> ";


        // Add all of the products
        for (ParticleType product : products){
            string += product.toString() + " + ";
        }

        // Remove the extra " + "
        string = string.substring(0, string.length() - 3);

        // Add the Q value
        string += String.format(" (%.2f MeV)", getQValue());

        // Newline
        string += "\n";

        // Return
        return string;
    }

    // Return cross section in units of cm^2
    public double getCrossSection(Particle p1, Particle p2){
        double centerMassEnergy = 0.0;

        centerMassEnergy += p1.getEnergy() * p2.getMass() / (p1.getMass() + p2.getMass());
        centerMassEnergy += p2.getEnergy() * p1.getMass() / (p1.getMass() + p2.getMass());

        try {
            return 1e-24*crossSection.value(centerMassEnergy);
        }catch (OutOfRangeException e){
            double[] knots = crossSection.getKnots();
            centerMassEnergy = FastMath.max(centerMassEnergy, knots[0]);
            centerMassEnergy = FastMath.min(centerMassEnergy, knots[knots.length-1]);
            return 1e-24*crossSection.value(centerMassEnergy);
        }
    }


    /**
     * Private convenience functions
     */

    public boolean equals(Object o){

        // Simple sanity checks
        if (o == this) return true;
        if (!(o instanceof NuclearReaction)) return false;

        // Verify that the reactants and products are the same length
        NuclearReaction reaction = (NuclearReaction) o;
        if (this.reactants.length != reaction.getReactants().length) return false;
        if (this.products.length  != reaction.getProducts().length)  return false;

        // Compare every reactant
        // TODO: Consider accounting for the fact that order doesn't matter?
        ParticleType[] reactants = reaction.getReactants();
        for (int i = 0; i < reactants.length; i++){
            if (!this.reactants[i].equals(reactants[i]))    return false;
        }

        // Compare every product
        // TODO: Consider accounting for the fact that order doesn't matter?
        ParticleType[] products = reaction.getProducts();
        for (int i = 0; i < products.length; i++){
            if (!this.products[i].equals(products[i]))    return false;
        }

        // If we're here, this is the same reaction
        return true;
    }


    // We'll use ZAID as our unique hash identifier
    public int hashCode(){
        int hashCode = 0;

        for (ParticleType reactant : reactants){
            hashCode += reactant.hashCode();
        }

        for (ParticleType product : products){
            hashCode += product.hashCode();
        }

        return hashCode;
    }

    private PolynomialSplineFunction generateCrossSectionFunction(File dataFile){

        ArrayList<Double> E = new ArrayList<>();
        ArrayList<Double> sig = new ArrayList<>();
        try{
            Scanner s = new Scanner(dataFile);

            while (s.hasNext()){
                E.add(s.nextDouble());
                sig.add(s.nextDouble());
            }

            double[] eArray = new double[E.size()];
            double[] sigArray = new double[sig.size()];
            for (int i = 0; i < E.size(); i++){
                eArray[i] = E.get(i);
                sigArray[i] = sig.get(i);
            }

            return new SplineInterpolator().interpolate(eArray, sigArray);
        }
        catch (IOException e){
            e.printStackTrace();
            return null;
        }
    }

    private Vector3D getCenterOfMassVelocity(Particle A, Particle B){
        Vector3D centerOfMassVelocity = A.getVelocity().scalarMultiply(A.getMass());                    // mA*vA
        centerOfMassVelocity = centerOfMassVelocity.add(B.getVelocity().scalarMultiply(B.getMass()));   // mA*vA + mB*vB
        centerOfMassVelocity = centerOfMassVelocity.scalarMultiply(1.0/(A.getMass() + B.getMass()));    // mA*vA + mB*vB / (mA + mB)
        return centerOfMassVelocity;
    }




}
