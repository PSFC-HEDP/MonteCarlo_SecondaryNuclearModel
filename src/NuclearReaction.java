import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.commons.math3.util.FastMath;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * Handles a nuclear reaction // A + B -> C + D
 */
public class NuclearReaction {

    private ParticleType reactantParticleA;
    private ParticleType reactantParticleB;

    private ParticleType productParticleC;
    private ParticleType productParticleD;

    private double energyReleased;      // Q value

    private PolynomialSplineFunction crossSection;

    public NuclearReaction(ParticleType A, ParticleType B, ParticleType C, ParticleType D, double Q, File crossSectionFile) {
        this.reactantParticleA = A;
        this.reactantParticleB = B;
        this.productParticleC = C;
        this.productParticleD = D;
        this.energyReleased = Q;

        generateCrossSectionFunction(crossSectionFile);
    }

    public Particle getProductParticleC(Particle A, Particle B){




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

    private double getCrossSection(Particle p1, Particle p2){
        double centerMassEnergy = 0.0;

        centerMassEnergy += p1.getE() * p2.getA() / (double) (p1.getA() + p2.getA());
        centerMassEnergy += p2.getE() * p1.getA() / (double) (p1.getA() + p2.getA());

        try {
            return crossSection.value(centerMassEnergy);
        }catch (OutOfRangeException e){
            double[] knots = crossSection.getKnots();
            centerMassEnergy = FastMath.max(centerMassEnergy, knots[0]);
            centerMassEnergy = FastMath.min(centerMassEnergy, knots[knots.length-1]);
            return crossSection.value(centerMassEnergy);
        }
    }


}
