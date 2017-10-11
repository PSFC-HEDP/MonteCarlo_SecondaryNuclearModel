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
 * Handles a nuclear reaction // A + B -> C + D
 */
public class NuclearReaction {

    private ParticleType reactantParticleTypeA;
    private ParticleType reactantParticleTypeB;

    private ParticleType productOfInterest;
    private ParticleType otherProduct;

    private double energyReleased;      // Q value

    private PolynomialSplineFunction crossSection;

    public static final NuclearReaction DTn   = new NuclearReaction(ParticleType.deuteron, ParticleType.triton, ParticleType.neutron, ParticleType.alpha, Constants.DT_N_ENERGY_RELEASE, Utils.DT_ENDF_XS_FILE);
    public static final NuclearReaction D3Hep = new NuclearReaction(ParticleType.deuteron, ParticleType.helium3, ParticleType.proton, ParticleType.alpha, Constants.D3HE_P_ENERGY_RELEASE, Utils.D3He_ENDF_XS_FILE);

    public NuclearReaction(ParticleType A, ParticleType B, ParticleType productOfInterest, ParticleType otherProduct, double Q, String crossSectionFile) {
        this.reactantParticleTypeA = A;
        this.reactantParticleTypeB = B;
        this.productOfInterest = productOfInterest;
        this.otherProduct = otherProduct;
        this.energyReleased = Q;

        generateCrossSectionFunction(new File(crossSectionFile));
    }

    public Particle getProductParticle(Particle A, Particle B, Vector3D direction){

        if (direction == null){
            return getProductParticle(A, B);
        }

        // Reduced mass of these two particles in amu
        double reducedMass = A.getMass()*B.getMass() / (A.getMass() + B.getMass());

        // Relative velocity of these two particles as a fraction of c
        Vector3D relativeVelocity = A.getVelocity().subtract(B.getVelocity());

        // Center of mass kinetic energy of this system in MeV
        double kineticEnergy = 0.5*Constants.MEV_PER_AMU*reducedMass*FastMath.pow(relativeVelocity.getNorm(),2);    // MeV

        // Center of mass velocity of this system in MeV
        Vector3D centerOfMassVelocity = getCenterOfMassVelocity(A, B);

        // Center of mass energy of the product of interest in MeV
        double productEnergy_CoM = (otherProduct.getMass() / (productOfInterest.getMass() + otherProduct.getMass())) * (energyReleased + kineticEnergy);

        // Center of mass speed of the product of interest as a fraction of c
        double productSpeed_CoM = FastMath.sqrt(2*productEnergy_CoM/productOfInterest.getMass()/Constants.MEV_PER_AMU);


        /**
         * How we handle it when we DO force the direction
         */

        // Angle between the center of mass velocity and the forced direction
        double angle = Vector3D.angle(centerOfMassVelocity, direction);

        // Use the law of cosines (based on vC = vCM + uC) to get the speed of the product as a fraction of c
        double productSpeed_Lab = FastMath.pow(centerOfMassVelocity.getNorm(), 2);       // vCM^2
        productSpeed_Lab *= FastMath.pow(FastMath.cos(angle),2) - 1;                     // vCM^2 * (cos^2(\theta) - 1)
        productSpeed_Lab += FastMath.pow(productSpeed_CoM,2);                    // uC^2 + vCM^2 * (cos^2(\theta) - 1)
        productSpeed_Lab = FastMath.sqrt(productSpeed_Lab);                                  // \sqrt(uC^2 + vCM^2 * (cos^2(\theta) - 1))
        productSpeed_Lab += centerOfMassVelocity.getNorm()*FastMath.cos(angle);             // vCM * cos(\theta) + \sqrt(uC^2 + vCM^2 * (cos^2(\theta) - 1))

        // Convert this to an energy in MeV
        double productEnergy_Lab = 0.5*Constants.MEV_PER_AMU*productOfInterest.getMass()*FastMath.pow(productSpeed_Lab,2);        // MeV


        // Build the  particle of type C born at the same position as A with the forced direction and our derived energy
        Particle productParticle = new Particle(productOfInterest, A.getPosition(), direction, productEnergy_Lab, A.getTime());

        return productParticle;
    }

    public Particle getProductParticle(Particle A, Particle B){

        // Reduced mass of these two particles in amu
        double reducedMass = A.getMass()*B.getMass() / (A.getMass() + B.getMass());

        // Relative velocity of these two particles as a fraction of c
        Vector3D relativeVelocity = A.getVelocity().subtract(B.getVelocity());

        // Center of mass kinetic energy of this system in MeV
        double kineticEnergy = 0.5*Constants.MEV_PER_AMU*reducedMass*FastMath.pow(relativeVelocity.getNorm(),2);    // MeV

        // Center of mass velocity of this system in MeV
        Vector3D centerOfMassVelocity = getCenterOfMassVelocity(A, B);

        // Center of mass energy of the product of interest in MeV
        double productEnergy_CoM = (otherProduct.getMass() / (productOfInterest.getMass() + otherProduct.getMass())) * (energyReleased + kineticEnergy);

        // Center of mass speed of the product of interest as a fraction of c
        double productSpeed_CoM = FastMath.sqrt(2*productEnergy_CoM/productOfInterest.getMass()/Constants.MEV_PER_AMU);


        /**
         * How we handle it if we didn't force the direction
         */

        Vector3D productVelocity_CoM = Utils.sampleRandomNormalizedVector().scalarMultiply(productSpeed_CoM);
        Vector3D productVelocity_Lab = centerOfMassVelocity.add(productVelocity_CoM);
        double productSpeed_Lab = productVelocity_Lab.getNorm();
        double productEnergy_Lab = 0.5*Constants.MEV_PER_AMU*productOfInterest.getMass()*FastMath.pow(productSpeed_Lab,2);        // MeV
        Particle productParticle = new Particle(productOfInterest, A.getPosition(), productVelocity_Lab.normalize(), productEnergy_Lab, A.getTime());


        return productParticle;
    }

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

    public ParticleType getReactantParticleTypeA() {
        return reactantParticleTypeA;
    }

    public ParticleType getReactantParticleTypeB() {
        return reactantParticleTypeB;
    }

    public ParticleType getProductOfInterest() {
        return productOfInterest;
    }

    public ParticleType getOtherProduct() {
        return otherProduct;
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
        if (o == this) return true;
        if (!(o instanceof NuclearReaction)) return false;

        NuclearReaction reaction = (NuclearReaction) o;
        return reactantParticleTypeA.equals(reaction.getReactantParticleTypeA())
                && reactantParticleTypeB.equals(reaction.getReactantParticleTypeB())
                && productOfInterest.equals(reaction.getProductOfInterest())
                && otherProduct.equals(reaction.getOtherProduct());
    }

    // We'll use ZAID as our unique hash identifier
    public int hashCode(){
        return reactantParticleTypeA.hashCode() +
                reactantParticleTypeB.hashCode() +
                productOfInterest.hashCode() +
                otherProduct.hashCode();
    }

    private void generateCrossSectionFunction(File dataFile){

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

            crossSection = new SplineInterpolator().interpolate(eArray, sigArray);
        }
        catch (IOException e){
            e.printStackTrace();
        }
    }

    private Vector3D getCenterOfMassVelocity(Particle A, Particle B){
        Vector3D centerOfMassVelocity = A.getVelocity().scalarMultiply(A.getMass());                    // mA*vA
        centerOfMassVelocity = centerOfMassVelocity.add(B.getVelocity().scalarMultiply(B.getMass()));   // mA*vA + mB*vB
        centerOfMassVelocity = centerOfMassVelocity.scalarMultiply(1.0/(A.getMass() + B.getMass()));    // mA*vA + mB*vB / (mA + mB)
        return centerOfMassVelocity;
    }

    public String toString(){
        return String.format("%s + %s -> %s + %s (%.2f MeV)",
                reactantParticleTypeA.toString(),
                reactantParticleTypeB.toString(),
                productOfInterest.toString(),
                otherProduct.toString(),
                energyReleased);
    }


}