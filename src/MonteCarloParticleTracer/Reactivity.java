package MonteCarloParticleTracer;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * Created by lahmann on 2016-06-16.
 */
public class Reactivity {

    private PolynomialSplineFunction value;


    public final static Reactivity DDn_Reactivity = ddn();
    public final static Reactivity DDp_Reactivity = ddp();
    public final static Reactivity DTn_Reactivity = dtn();
    public final static Reactivity D3Hep_Reactivity = d3Hep();
    public final static Reactivity He3He3_Reactivity = he3he3();


    private static Reactivity ddn(){
        return new Reactivity(DataFiles.DDn_Reactivity_File);
    }
    private static Reactivity ddp(){
        return new Reactivity(DataFiles.DDp_Reactivity_File);
    }
    private static Reactivity dtn() { return new Reactivity(DataFiles.DTn_Reactivity_File); }
    private static Reactivity d3Hep(){
        return new Reactivity(DataFiles.D3Hep_Reactivity_File);
    }
    private static Reactivity he3he3() {
        return new Reactivity(DataFiles.He3He3_Reactivity_File);
    }



    public Reactivity(File dataFile){

        ArrayList<Double> T = new ArrayList<>();
        ArrayList<Double> sigV = new ArrayList<>();
        try{
            Scanner s = new Scanner(dataFile);

            while (s.hasNext()){
                T.add(s.nextDouble());
                sigV.add(s.nextDouble());
            }

            double[] eArray = new double[T.size()];
            double[] sigArray = new double[sigV.size()];
            for (int i = 0; i < T.size(); i++){
                eArray[i] = T.get(i);
                sigArray[i] = sigV.get(i);
            }

            value = new LinearInterpolator().interpolate(
                    eArray, sigArray);
        }
        catch (IOException e){
            e.printStackTrace();
        }
    }


    double evaluate(double temperature){
        return value.value(temperature);
    }
}
