import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * Created by lahmann on 2016-06-16.
 */
public class CrossSection {

    private PolynomialSplineFunction function;

    public CrossSection(String dataFilename){
        this(new File(dataFilename));
    }

    public CrossSection(File dataFile){

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

            function = new SplineInterpolator().interpolate(
                    eArray, sigArray);
        }
        catch (IOException e){
            e.printStackTrace();
        }
    }

    // TODO: This is currently assuming the limit v1 >> v2
    public double evaluate(Particle fastParticle, Particle slowParticle){
        double centerMassEnergy = (double) slowParticle.getA() / (fastParticle.getA() + slowParticle.getA());
        centerMassEnergy *= fastParticle.getE();
        return function.value(centerMassEnergy);
    }

    public static CrossSection dt(){
        return new CrossSection(Utils.DT_ENDF_XS_FILE);
    }

    public static CrossSection d3He(){
        return new CrossSection(Utils.D3He_ENDF_XS_FILE);
    }
}
